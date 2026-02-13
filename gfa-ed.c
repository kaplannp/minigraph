#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include "gfa-priv.h"
#include "kalloc.h"
#include "ksort.h"
#include "khashl.h" // make it compatible with kalloc
#include "kdq.h"
#include "kvec-km.h"

int gfa_ed_dbg = 0;

static FILE *gwf_ed_stats_fp(void) {
	static FILE *fp;
	if (!fp) fp = fopen("gwf-ed-stats.log", "a");
	return fp;
}

/***************
 * Preparation *
 ***************/

void gfa_edopt_init(gfa_edopt_t *opt)
{
	memset(opt, 0, sizeof(gfa_edopt_t));
	opt->bw_dyn = opt->max_lag = opt->s_term = -1;
	opt->max_chk = 1000;
}

gfa_edseq_t *gfa_edseq_init(const gfa_t *g)
{
	uint32_t i, n_vtx = gfa_n_vtx(g);
	gfa_edseq_t *es;
	GFA_MALLOC(es, n_vtx);
	for (i = 0; i < g->n_seg; ++i) {
		const gfa_seg_t *s = &g->seg[i];
		char *t;
		int32_t j;
		GFA_MALLOC(t, s->len + 1);
		for (j = 0; j < s->len; ++j)
			t[s->len - j - 1] = gfa_comp_table[(uint8_t)s->seq[j]];
		t[s->len] = 0;
		es[i<<1].seq = (char*)s->seq;
		es[i<<1|1].seq = t;
		es[i<<1].len = es[i<<1|1].len = s->len;
	}
	return es;
}

void gfa_edseq_destroy(int32_t n_seg, gfa_edseq_t *es)
{
	int32_t i;
	for (i = 0; i < n_seg; ++i)
		free((char*)es[i<<1|1].seq);
	free(es);
}

/*****************
 * Edit distance *
 *****************/

#define GWF_DIAG_SHIFT 0x40000000

static inline uint64_t gwf_gen_vd(uint32_t v, int32_t d)
{
	return (uint64_t)v<<32 | (GWF_DIAG_SHIFT + d);
}

/*
 * Diagonal interval
 */
typedef struct {
	uint64_t vd0, vd1;
} gwf_intv_t;

typedef kvec_t(gwf_intv_t) gwf_intv_v;

#define intvd_key(x) ((x).vd0)
KRADIX_SORT_INIT(gwf_intv, gwf_intv_t, intvd_key, 8)

static int gwf_intv_is_sorted(int32_t n_a, const gwf_intv_t *a)
{
	int32_t i;
	for (i = 1; i < n_a; ++i)
		if (a[i-1].vd0 > a[i].vd0) break;
	return (i == n_a);
}

// merge overlapping intervals; input must be sorted
static size_t gwf_intv_merge_adj(size_t n, gwf_intv_t *a)
{
	size_t i, k;
	uint64_t st, en;
	if (n == 0) return 0;
	st = a[0].vd0, en = a[0].vd1;
	for (i = 1, k = 0; i < n; ++i) {
		if (a[i].vd0 > en) {
			a[k].vd0 = st, a[k++].vd1 = en;
			st = a[i].vd0, en = a[i].vd1;
		} else en = en > a[i].vd1? en : a[i].vd1;
	}
	a[k].vd0 = st, a[k++].vd1 = en;
	return k;
}

// merge two sorted interval lists
static size_t gwf_intv_merge2(gwf_intv_t *a, size_t n_b, const gwf_intv_t *b, size_t n_c, const gwf_intv_t *c)
{
	size_t i = 0, j = 0, k = 0;
	while (i < n_b && j < n_c) {
		if (b[i].vd0 <= c[j].vd0)
			a[k++] = b[i++];
		else a[k++] = c[j++];
	}
	while (i < n_b) a[k++] = b[i++];
	while (j < n_c) a[k++] = c[j++];
	return gwf_intv_merge_adj(k, a);
}

/*
 * Diagonal
 * Observe that there are no arrays. This implies that it must be a single
 * diagonal not a wavefront
 */
typedef struct { // a diagonal
	uint64_t vd; // higher 32 bits: vertex ID; lower 32 bits: diagonal+0x4000000
	int32_t k; //this is how much we've pushed down diagonal I think (d in paper)
	int32_t len;
	uint32_t xo; // higher 31 bits: anti diagonal; lower 1 bit: out-of-order or not
	int32_t t; //By process of elimination this is score? Or they have no score
} gwf_diag_t;

typedef kvec_t(gwf_diag_t) gwf_diag_v;

#define ed_key(x) ((x).vd)
KRADIX_SORT_INIT(gwf_ed, gwf_diag_t, ed_key, 8)

KDQ_INIT(gwf_diag_t)

// push (v,d,k) to the end of the queue
static inline void gwf_diag_push(void *km, gwf_diag_v *a, uint32_t v, int32_t d, int32_t k, uint32_t x, uint32_t ooo, int32_t t)
{
	gwf_diag_t *p;
	kv_pushp(gwf_diag_t, km, *a, &p);
	p->vd = gwf_gen_vd(v, d), p->k = k, p->xo = x<<1|ooo, p->t = t;
}

// determine the wavefront on diagonal (v,d)
//zkn ok, so this d is what we would usually call k from the paper. It is the
//index of the diagonal.  In contrast k is the distance which we would usually
//call d in the paper.
static inline int32_t gwf_diag_update(gwf_diag_t *p, uint32_t v, int32_t d, int32_t k, uint32_t x, uint32_t ooo, int32_t t)
{
	uint64_t vd = gwf_gen_vd(v, d);
	if (p->vd == vd) {
		p->xo = p->k > k? p->xo : x<<1|ooo;
		p->t  = p->k > k? p->t : t;
		p->k  = p->k > k? p->k : k; 
		return 0;
	}
	return 1;
}

static int gwf_diag_is_sorted(int32_t n_a, const gwf_diag_t *a)
{
	int32_t i;
	for (i = 1; i < n_a; ++i)
		if (a[i-1].vd > a[i].vd) break;
	return (i == n_a);
}

// sort a[]. This uses the gwf_diag_t::ooo field to speed up sorting.
static void gwf_diag_sort(int32_t n_a, gwf_diag_t *a, void *km, gwf_diag_v *ooo)
{
	int32_t i, j, k, n_b, n_c;
	gwf_diag_t *b, *c;

	kv_resize(gwf_diag_t, km, *ooo, n_a);
	for (i = 0, n_c = 0; i < n_a; ++i)
		if (a[i].xo&1) ++n_c;
	n_b = n_a - n_c;
	b = ooo->a, c = b + n_b;
	for (i = j = k = 0; i < n_a; ++i) {
		if (a[i].xo&1) c[k++] = a[i];
		else b[j++] = a[i];
	}
	radix_sort_gwf_ed(c, c + n_c);
	for (k = 0; k < n_c; ++k) c[k].xo &= 0xfffffffeU;

	i = j = k = 0;
	while (i < n_b && j < n_c) {
		if (b[i].vd <= c[j].vd)
			a[k++] = b[i++];
		else a[k++] = c[j++];
	}
	while (i < n_b) a[k++] = b[i++];
	while (j < n_c) a[k++] = c[j++];
}

// remove diagonals not on the wavefront
static int32_t gwf_diag_dedup(int32_t n_a, gwf_diag_t *a, void *km, gwf_diag_v *ooo)
{
	int32_t i, n, st;
	if (!gwf_diag_is_sorted(n_a, a))
		gwf_diag_sort(n_a, a, km, ooo);
	for (i = 1, st = 0, n = 0; i <= n_a; ++i) {
		if (i == n_a || a[i].vd != a[st].vd) {
			int32_t j, max_j = st;
			if (st + 1 < i)
				for (j = st + 1; j < i; ++j) // choose the far end (i.e. the wavefront)
					if (a[max_j].k < a[j].k) max_j = j;
			a[n++] = a[max_j];
			st = i;
		}
	}
	return n;
}

// use forbidden bands to remove diagonals not on the wavefront
static int32_t gwf_mixed_dedup(int32_t n_a, gwf_diag_t *a, int32_t n_b, gwf_intv_t *b)
{
	int32_t i = 0, j = 0, k = 0;
	while (i < n_a && j < n_b) {
		if (a[i].vd >= b[j].vd0 && a[i].vd < b[j].vd1) ++i;
		else if (a[i].vd >= b[j].vd1) ++j;
		else a[k++] = a[i++];
	}
	while (i < n_a) a[k++] = a[i++];
	return k;
}

/*
 * Traceback stack
 */
KHASHL_MAP_INIT(KH_LOCAL, gwf_map64_t, gwf_map64, uint64_t, int32_t, kh_hash_uint64, kh_eq_generic)

typedef struct {
	int32_t v;
	int32_t pre;
} gwf_trace_t;

typedef kvec_t(gwf_trace_t) gwf_trace_v;

static int32_t gwf_trace_push(void *km, gwf_trace_v *a, int32_t v, int32_t pre, gwf_map64_t *h)
{
	uint64_t key = (uint64_t)v << 32 | (uint32_t)pre;
	khint_t k;
	int absent;
	k = gwf_map64_put(h, key, &absent);
	if (absent) {
		gwf_trace_t *p;
		kv_pushp(gwf_trace_t, km, *a, &p);
		p->v = v, p->pre = pre;
		kh_val(h, k) = a->n - 1;
		return a->n - 1;
	}
	return kh_val(h, k);
}

/*
 * Core GWFA routine
 */
KHASHL_INIT(KH_LOCAL, gwf_set64_t, gwf_set64, uint64_t, kh_hash_dummy, kh_eq_generic)

typedef struct {
	void *km;
	gwf_set64_t *ha; // hash table for adjacency
	gwf_map64_t *ht; // hash table for traceback
	gwf_intv_v intv;
	gwf_intv_v tmp, swap;
	gwf_diag_v ooo;
	gwf_trace_v t;
} gwf_edbuf_t;

// remove diagonals not on the wavefront
static int32_t gwf_dedup(gwf_edbuf_t *buf, int32_t n_a, gwf_diag_t *a)
{
	if (buf->intv.n + buf->tmp.n > 0) {
		if (!gwf_intv_is_sorted(buf->tmp.n, buf->tmp.a))
			radix_sort_gwf_intv(buf->tmp.a, buf->tmp.a + buf->tmp.n);
		kv_copy(gwf_intv_t, buf->km, buf->swap, buf->intv);
		kv_resize(gwf_intv_t, buf->km, buf->intv, buf->intv.n + buf->tmp.n);
		buf->intv.n = gwf_intv_merge2(buf->intv.a, buf->swap.n, buf->swap.a, buf->tmp.n, buf->tmp.a);
	}
	n_a = gwf_diag_dedup(n_a, a, buf->km, &buf->ooo);
	if (buf->intv.n > 0)
		n_a = gwf_mixed_dedup(n_a, a, buf->intv.n, buf->intv.a);
	return n_a;
}

// remove diagonals that lag far behind the furthest wavefront
static int32_t gwf_prune(int32_t n_a, gwf_diag_t *a, uint32_t max_lag, int32_t bw_dyn)
{
	int32_t i, j, iq, dq, max_i = -1;
	uint32_t max_x = 0;
	gwf_diag_t *q;
	for (i = 0; i < n_a; ++i)
		if (a[i].xo>>1 > max_x)
			max_x = a[i].xo>>1, max_i = i;
	q = &a[max_i];
	iq = (int32_t)q->vd - GWF_DIAG_SHIFT + q->k;
	dq = (int32_t)(q->xo>>1) - iq - iq;
	for (i = j = 0; i < n_a; ++i) {
		gwf_diag_t *p = &a[i];
		int32_t ip = (int32_t)p->vd - GWF_DIAG_SHIFT + p->k;
		int32_t dp = (int32_t)(p->xo>>1) - ip - ip;
		int32_t w = dp > dq? dp - dq : dq - dp;
		if (bw_dyn >= 0 && w > bw_dyn) continue;
		if ((p->xo>>1) + max_lag < max_x) continue;
		a[j++] = *p;
	}
	return j;
}

//zkn this is the function that pushes the wavefront along the diagonal with
//matches. We do examine the query here
// reach the wavefront
static inline int32_t gwf_extend1(int32_t d, int32_t k, int32_t vl, const char *ts, int32_t ql, const char *qs)
{
	int32_t max_k = (ql - d < vl? ql - d : vl) - 1;
	const char *ts_ = ts + 1, *qs_ = qs + d + 1;
#if 0
	// int32_t i = k + d; while (k + 1 < vl && i + 1 < ql && ts[k+1] == q[i+1]) ++k, ++i;
	while (k < max_k && *(ts_ + k) == *(qs_ + k))
		++k;
#else
	uint64_t cmp = 0;
	while (k + 7 < max_k) {
		uint64_t x = *(uint64_t*)(ts_ + k); // warning: unaligned memory access
		uint64_t y = *(uint64_t*)(qs_ + k);
		cmp = x ^ y;
		if (cmp == 0) k += 8;
		else break;
	}
	if (cmp)
		k += __builtin_ctzl(cmp) >> 3; // on x86, this is done via the BSR instruction: https://www.felixcloutier.com/x86/bsr
	else if (k + 7 >= max_k)
		while (k < max_k && *(ts_ + k) == *(qs_ + k)) // use this for generic CPUs. It is slightly faster than the unoptimized version
			++k;
#endif
	return k;
}

// This is essentially Landau-Vishkin for linear sequences. The function speeds up alignment to long vertices. Not really necessary.
static void gwf_ed_extend_batch(void *km, const gfa_t *g, const gfa_edseq_t *es, int32_t ql, const char *q, int32_t n, gwf_diag_t *a, gwf_diag_v *B,
								kdq_t(gwf_diag_t) *A, gwf_intv_v *tmp_intv, gfa_edrst_t *r)
{
	int32_t j, m;
	int32_t v = a->vd>>32;
	int32_t vl = es[v].len;
	const char *ts = es[v].seq;
	gwf_diag_t *b;

	// wfa_extend
	for (j = 0; j < n; ++j) {
		int32_t k;
		k = gwf_extend1((int32_t)a[j].vd - GWF_DIAG_SHIFT, a[j].k, vl, ts, ql, q);
		a[j].len = k - a[j].k;
		a[j].xo += a[j].len << 2;
		a[j].k = k;
	}

	// wfa_next
	kv_resize(gwf_diag_t, km, *B, B->n + n + 2);
	b = &B->a[B->n];
	b[0].vd = a[0].vd - 1;
	b[0].xo = a[0].xo + 2; // 2 == 1<<1
	b[0].k = a[0].k + 1;
	b[0].t = a[0].t;
	b[1].vd = a[0].vd;
	b[1].xo =  n == 1 || a[0].k > a[1].k? a[0].xo + 4 : a[1].xo + 2;
	b[1].t  =  n == 1 || a[0].k > a[1].k? a[0].t : a[1].t;
	b[1].k  = (n == 1 || a[0].k > a[1].k? a[0].k : a[1].k) + 1;
	for (j = 1; j < n - 1; ++j) {
		uint32_t x = a[j-1].xo + 2;
		int32_t k = a[j-1].k, t = a[j-1].t;
		x = k > a[j].k + 1? x : a[j].xo + 4;
		t = k > a[j].k + 1? t : a[j].t;
		k = k > a[j].k + 1? k : a[j].k + 1;
		x = k > a[j+1].k + 1? x : a[j+1].xo + 2;
		t = k > a[j+1].k + 1? t : a[j+1].t;
		k = k > a[j+1].k + 1? k : a[j+1].k + 1;
		b[j+1].vd = a[j].vd, b[j+1].k = k, b[j+1].xo = x, b[j+1].t = t;
	}
	if (n >= 2) {
		b[n].vd = a[n-1].vd;
		b[n].xo = a[n-2].k > a[n-1].k + 1? a[n-2].xo + 2 : a[n-1].xo + 4;
		b[n].t  = a[n-2].k > a[n-1].k + 1? a[n-2].t : a[n-1].t;
		b[n].k  = a[n-2].k > a[n-1].k + 1? a[n-2].k : a[n-1].k + 1;
	}
	b[n+1].vd = a[n-1].vd + 1;
	b[n+1].xo = a[n-1].xo + 2;
	b[n+1].t  = a[n-1].t;
	b[n+1].k  = a[n-1].k;

	// drop out-of-bound cells
	//if (a[n-1].k == vl - 1) b[n+1].k = vl; // insertion to the end of a vertex is handled elsewhere. FIXME: this line leads to wrong result for MHC-57 and MHC-HG002.2
	for (j = 0; j < n; ++j) {
		gwf_diag_t *p = &a[j];
		if (p->k == vl - 1 || (int32_t)p->vd - GWF_DIAG_SHIFT + p->k == ql - 1)
			p->xo |= 1, *kdq_pushp(gwf_diag_t, A) = *p;
	}
	for (j = 0, m = 0; j < n + 2; ++j) {
		gwf_diag_t *p = &b[j];
		int32_t d = (int32_t)p->vd - GWF_DIAG_SHIFT;
		if (d + p->k < ql && p->k < vl) {
			b[m++] = *p;
		} else if (p->k == vl) {
			gwf_intv_t *q;
			kv_pushp(gwf_intv_t, km, *tmp_intv, &q);
			q->vd0 = gwf_gen_vd(v, d), q->vd1 = q->vd0 + 1;
		}
	}
	B->n += m;
}

//This function demonstrates, 1) that we are not doing affine gap, and 2) that
//we are doing base level alignment
// wfa_extend and wfa_next combined
static gwf_diag_t *gwf_ed_extend(gwf_edbuf_t *buf, const gfa_edopt_t *opt, const gfa_t *g, const gfa_edseq_t *es, int32_t s, int32_t ql, const char *q,
								 uint32_t v1, int32_t off1, int32_t *end_tb, int32_t *n_a_, gwf_diag_t *a, gfa_edrst_t *r)
{
	int32_t i, x, n = *n_a_, do_dedup = 1;
	kdq_t(gwf_diag_t) *A; //A is the queue containing the verticies to look at!
	gwf_diag_v B = {0,0,0};
	gwf_diag_t *b;
	FILE *stats_fp = gwf_ed_stats_fp();

	if (stats_fp) {
		fprintf(stats_fp, "SCORE: %d\n", s);
		fprintf(stats_fp, "n: %d\n", n);
		fflush(stats_fp);
	}

	r->end_v = (uint32_t)-1;
	r->end_off = *end_tb = -1;
	buf->tmp.n = 0;
	gwf_set64_clear(buf->ha); // hash table $h to avoid visiting a vertex twice
	for (i = 0, x = 1; i < 32; ++i, x <<= 1)
		if (x >= n) break;
	if (i < 4) i = 4;
	A = kdq_init2(gwf_diag_t, buf->km, i); // $A is a queue
	kv_resize(gwf_diag_t, buf->km, B, n * 2);
#if 0 // unoptimized version without calling gwf_ed_extend_batch() at all. The final result will be the same.
	A->count = n;
	memcpy(A->a, a, n * sizeof(*a));
#else // optimized for long vertices.
	for (x = 0, i = 1; i <= n; ++i) {
		if (i == n || a[i].vd != a[i-1].vd + 1) {
			gwf_ed_extend_batch(buf->km, g, es, ql, q, i - x, &a[x], &B, A, &buf->tmp, r);
			if (stats_fp) {
				//fprintf(stats_fp, "EXTEND_BATCH: %u COUNT: %d\n", (uint32_t)(a[x].vd >> 32), i - x);
				fflush(stats_fp);
			}
			x = i;
		}
	}
	if (kdq_size(A) == 0) do_dedup = 0;
#endif
	kfree(buf->km, a); // $a is not used as it has been copied to $A

	fprintf(stats_fp, "SIZE A: %d\n", (int)kdq_size(A));
	while (kdq_size(A)) {
		gwf_diag_t t; 
		uint32_t v, x0;
		int32_t ooo, d, k, i, vl;

		t = *kdq_shift(gwf_diag_t, A);
		ooo = t.xo&1, v = t.vd >> 32; // vertex
		d = (int32_t)t.vd - GWF_DIAG_SHIFT; // diagonal
		k = t.k; // wavefront position on the vertex
		vl = es[v].len; // $vl is the vertex length
		k = gwf_extend1(d, k, vl, es[v].seq, ql, q); //push along diagonal far as you can
		i = k + d; // query position
		x0 = (t.xo >> 1) + ((k - t.k) << 1); // current anti diagonal

		if (k + 1 < vl && i + 1 < ql) { // the most common case: the wavefront is in the middle
			int32_t push1 = 1, push2 = 1;
      //So somehow, B.a appears to hold the wavefront
      // I believe these first two are the cases of the switch statement. They
      // are the initial updates from previous scores that is. (and there isn't
      // any affine gap)
			if (B.n >= 2) push1 = gwf_diag_update(&B.a[B.n - 2], v, d-1, k+1, x0 + 1, ooo, t.t);
			if (B.n >= 1) push2 = gwf_diag_update(&B.a[B.n - 1], v, d,   k+1, x0 + 2, ooo, t.t);
      //Now we add it to the queue. When did we extend it?
			if (push1)          gwf_diag_push(buf->km, &B, v, d-1, k+1, x0 + 1, 1, t.t);
			if (push2 || push1) gwf_diag_push(buf->km, &B, v, d,   k+1, x0 + 2, 1, t.t);
			gwf_diag_push(buf->km, &B, v, d+1, k, x0 + 1, ooo, t.t);
    //This else is the branchy weirdness part where you keep pushing along edges
		} else if (i + 1 < ql) { // k + 1 == g->len[v]; reaching the end of the vertex but not the end of query
			int32_t nv = gfa_arc_n(g, v), j, n_ext = 0, tw = -1;
			gfa_arc_t *av = gfa_arc_a(g, v);
			gwf_intv_t *p;
			kv_pushp(gwf_intv_t, buf->km, buf->tmp, &p);
			p->vd0 = gwf_gen_vd(v, d), p->vd1 = p->vd0 + 1;
			if (opt->traceback) tw = gwf_trace_push(buf->km, &buf->t, v, t.t, buf->ht);
			for (j = 0; j < nv; ++j) { // traverse $v's neighbors
				uint32_t w = av[j].w; // $w is next to $v
				int32_t ol = av[j].ow;
				int absent;
				gwf_set64_put(buf->ha, (uint64_t)w<<32 | (i + 1), &absent); // test if ($w,$i) has been visited
				if (q[i + 1] == es[w].seq[ol]) { // can be extended to the next vertex without a mismatch
					++n_ext;
					if (absent) {
						gwf_diag_t *p;
						p = kdq_pushp(gwf_diag_t, A);
						p->vd = gwf_gen_vd(w, i + 1 - ol), p->k = ol, p->xo = (x0+2)<<1 | 1, p->t = tw;
						//if (stats_fp) {
						//	fprintf(stats_fp, "TRANSITION_A_V: %u A_SIZE: %d\n", w, (int)kdq_size(A));
						//	fflush(stats_fp);
						//}
					}
				} else if (absent) {
					gwf_diag_push(buf->km, &B, w, i - ol,     ol, x0 + 1, 1, tw);
					gwf_diag_push(buf->km, &B, w, i + 1 - ol, ol, x0 + 2, 1, tw);
				}
			}
			if (nv == 0 || n_ext != nv) // add an insertion to the target; this *might* cause a duplicate in corner cases
				gwf_diag_push(buf->km, &B, v, d+1, k, x0 + 1, 1, t.t);
		} else if (v1 == (uint32_t)-1 || (v == v1 && k == off1)) { // i + 1 == ql
			r->end_v = v, r->end_off = k, r->wlen = x0 - i - 1, *end_tb = t.t, *n_a_ = 0;
			if (stats_fp) {
				fprintf(stats_fp, "TARGET_REACHED: %u\n", v);
				fflush(stats_fp);
			}
			kdq_destroy(gwf_diag_t, A);
			kfree(buf->km, B.a);
			return 0;
		} else if (k + 1 < vl) { // i + 1 == ql; reaching the end of the query but not the end of the vertex
			gwf_diag_push(buf->km, &B, v, d-1, k+1, x0 + 1, ooo, t.t); // add an deletion; this *might* case a duplicate in corner cases
		} else if (v != v1) { // i + 1 == ql && k + 1 == g->len[v]; not reaching the last vertex $v1
			int32_t nv = gfa_arc_n(g, v), j, tw = -1;
			const gfa_arc_t *av = gfa_arc_a(g, v);
			if (opt->traceback) tw = gwf_trace_push(buf->km, &buf->t, v, t.t, buf->ht);
			for (j = 0; j < nv; ++j)
				gwf_diag_push(buf->km, &B, av[j].w, i - av[j].ow, av[j].ow, x0 + 1, 1, tw); // deleting the first base on the next vertex
		} else { // may come here when k>off1 (due to banding); do nothing in this case
		}
	}

	kdq_destroy(gwf_diag_t, A);
	*n_a_ = n = B.n, b = B.a;

	if (stats_fp) {
		fprintf(stats_fp, "N_A_NEXT: %d\n", n);
		fflush(stats_fp);
	}

	if (do_dedup) *n_a_ = n = gwf_dedup(buf, n, b);
	if (opt->max_lag > 0 && n > opt->max_chk && ((s+1)&0xf) == 0)
		*n_a_ = n = gwf_prune(n, b, opt->max_lag, opt->bw_dyn);
	if (stats_fp) {
		int64_t n_inactive = 0;
		size_t ii;
		for (ii = 0; ii < buf->intv.n; ++ii)
			n_inactive += buf->intv.a[ii].vd1 - buf->intv.a[ii].vd0;
		fprintf(stats_fp, "N_INACTIVE: %lld\n", (long long)n_inactive);
		fprintf(stats_fp, "INTV_BYTES: %lld\n",
			(long long)(buf->intv.n * sizeof(gwf_intv_t)));
		fflush(stats_fp);
	}
	return b;
}

static void gwf_traceback(gwf_edbuf_t *buf, int32_t end_v, int32_t end_tb, gfa_edrst_t *path)
{
	int32_t i = end_tb, n = 1;
	while (i >= 0 && buf->t.a[i].v >= 0)
		++n, i = buf->t.a[i].pre;
	KMALLOC(buf->km, path->v, n);
	i = end_tb, n = 0;
	path->v[n++] = end_v;
	while (i >= 0 && buf->t.a[i].v >= 0)
		path->v[n++] = buf->t.a[i].v, i = buf->t.a[i].pre;
	path->nv = n;
	for (i = 0; i < path->nv>>1; ++i)
		n = path->v[i], path->v[i] = path->v[path->nv - 1 - i], path->v[path->nv - 1 - i] = n;
}

static void gwf_ed_print_diag(const gfa_t *g, size_t n, gwf_diag_t *a) // for debugging only
{
	size_t i;
	for (i = 0; i < n; ++i) {
		int32_t d = (int32_t)a[i].vd - GWF_DIAG_SHIFT;
		printf("Z\t%d\t%s\t%d\t%d\t%d\n", d + a[i].k, g->seg[(a[i].vd>>32)>>1].name, d, a[i].k, a[i].xo>>1);
	}
}

static void gwf_ed_print_intv(size_t n, gwf_intv_t *a) // for debugging only
{
	size_t i;
	for (i = 0; i < n; ++i)
		printf("Z\t%d\t%d\t%d\n", (int32_t)(a[i].vd0>>32), (int32_t)a[i].vd0 - GWF_DIAG_SHIFT, (int32_t)a[i].vd1 - GWF_DIAG_SHIFT);
}

/*
 * BFS to find all vertices on any valid path from v0 to v1.
 *
 * A "valid path" is one whose total sequence length <= budget
 * (where budget = query gap length, i.e. z->ql).
 *
 * Algorithm: two-pass weighted BFS + intersection.
 *
 *   Pass 1 (forward): BFS from v0 along outgoing arcs.
 *     Computes fwd[v] = min total sequence length on any
 *     path from v0 to v (including both endpoints).
 *     NOTE* BFS, but we're doing unweighted. so we still need to check when we revist a node it 
 *     could be shorter
 *
 *   Pass 2 (backward): BFS from v1 along incoming arcs.
 *     Computes bwd[v] = min total sequence length on any
 *     path from v to v1 (including both endpoints).
 *     Incoming arcs are found via the complement convention:
 *       arc v->w  <=>  complement arc w^1->v^1
 *     So outgoing arcs from v^1 give us the predecessors
 *     of v (each target w maps to predecessor w^1).
 *
 *   Intersection: vertex v is on a valid path iff
 *     fwd[v] + bwd[v] - es[v].len <= budget
 *     (subtract es[v].len because v's length is counted in
 *     both fwd and bwd distances).
 *
 * Both passes use relaxation-BFS: a vertex is re-enqueued
 * whenever a shorter distance is found. This is correct for
 * positive weights (all sequence lengths > 0), similar to
 * Bellman-Ford but with a FIFO queue.
 *
 * @param km      thread-local memory allocator
 * @param g       the graph (for arc traversal)
 * @param es      per-vertex sequences (for lengths)
 * @param v0      source vertex (start of GWFA alignment)
 * @param v1      sink vertex (target of GWFA alignment)
 * @param budget  max total sequence length (= query length)
 * @param result  output: set of vertices on valid paths
 */
static void gwf_bfs_reachable(void *km, const gfa_t *g,
	const gfa_edseq_t *es, uint32_t v0, uint32_t v1,
	int32_t budget, gwf_set64_t *result)
{
	gwf_map64_t *fwd, *bwd;      // distance maps: vertex -> min dist
	typedef kvec_t(uint32_t) uint32_v;
	uint32_v queue = {0, 0, 0};  // BFS queue (reused for both passes)
	int32_t head, j;              // head = queue read pointer
	khint_t k;                    // hash table bucket index
	int absent;                   // set by put: 1 if key was new

	/*
	 * PASS 1: Forward BFS from v0
	 *
	 * fwd maps vertex_id (uint64_t) -> min distance (int32_t).
	 * Distance = sum of sequence lengths along the path,
	 * including the source vertex v0.
	 *
	 * We use gwf_map64_t (KHASHL_MAP_INIT with uint64_t keys
	 * and int32_t values). Access pattern:
	 *   gwf_map64_put(h, key, &absent) -> bucket index
	 *   gwf_map64_get(h, key) -> bucket index (kh_end if missing)
	 *   kh_val(h, bucket) -> the int32_t value
	 */
	fwd = gwf_map64_init2(km);
	// Seed: distance to v0 = its own sequence length
	k = gwf_map64_put(fwd, (uint64_t)v0, &absent);
	kh_val(fwd, k) = es[v0].len;
	kv_push(uint32_t, km, queue, v0);

	// Process queue. Note: queue.n grows as we enqueue new
	// vertices, so this loop processes all reachable vertices.
	for (head = 0; head < (int32_t)queue.n; ++head) {
		uint32_t v = queue.a[head];
		// Look up current min distance to v
		khint_t kv = gwf_map64_get(fwd, (uint64_t)v);
		int32_t dv = kh_val(fwd, kv);
		// Get outgoing arcs from v: gfa_arc_a returns the arc
		// array, gfa_arc_n returns the count
		int32_t nv = gfa_arc_n(g, v);
		const gfa_arc_t *av = gfa_arc_a(g, v);
		for (j = 0; j < nv; ++j) {
			uint32_t w = av[j].w;       // neighbor vertex
			// Candidate distance to w = dist to v + w's length
			int32_t dw = dv + es[w].len;
			khint_t kw;
			if (dw > budget) continue;  // over budget, prune
			// Insert or find w in the distance map.
			// If absent=1, w is new. If absent=0, w exists.
			kw = gwf_map64_put(fwd, (uint64_t)w, &absent);
			// Update if this is a new vertex or we found a
			// shorter path than previously recorded.
			if (absent || kh_val(fwd, kw) > dw) {
				kh_val(fwd, kw) = dw;
				// Re-enqueue for further exploration with
				// the shorter distance (relaxation).
				kv_push(uint32_t, km, queue, w);
			}
		}
	}

	/*
	 * PASS 2: Backward BFS from v1 (reverse edges)
	 *
	 * bwd maps vertex_id -> min distance from that vertex to v1.
	 *
	 * To traverse incoming edges (predecessors), we use the
	 * complement arc convention from gfa.h:
	 *   If arc (a -> v) exists, its complement (v^1 -> a^1) exists.
	 * So: outgoing arcs from v^1 have targets {w0, w1, ...},
	 *   and the actual predecessors of v are {w0^1, w1^1, ...}.
	 * (^1 flips the orientation bit: v^1 = v XOR 1)
	 */
	bwd = gwf_map64_init2(km);
	// Seed: distance from v1 to itself = its own sequence length
	k = gwf_map64_put(bwd, (uint64_t)v1, &absent);
	kh_val(bwd, k) = es[v1].len;
	// Reuse the queue array, reset length
	queue.n = 0;
	kv_push(uint32_t, km, queue, v1);

	for (head = 0; head < (int32_t)queue.n; ++head) {
		uint32_t v = queue.a[head];
		khint_t kv = gwf_map64_get(bwd, (uint64_t)v);
		int32_t dv = kh_val(bwd, kv);
		// Get predecessors of v via complement arcs:
		// outgoing arcs from v^1, each target w -> pred = w^1
		int32_t nv = gfa_arc_n(g, v ^ 1);
		const gfa_arc_t *av = gfa_arc_a(g, v ^ 1);
		for (j = 0; j < nv; ++j) {
			uint32_t pred = av[j].w ^ 1; // actual predecessor
			// Distance from pred to v1 = dist(v->v1) + pred's len
			int32_t dp = dv + es[pred].len;
			khint_t kp;
			if (dp > budget) continue;  // over budget, prune
			kp = gwf_map64_put(bwd, (uint64_t)pred, &absent);
			if (absent || kh_val(bwd, kp) > dp) {
				kh_val(bwd, kp) = dp;
				kv_push(uint32_t, km, queue, pred);
			}
		}
	}

	/*
	 * INTERSECTION: find vertices on valid paths.
	 *
	 * A vertex v is on some path v0->...->v->...->v1 with total
	 * sequence length <= budget iff:
	 *   fwd[v] + bwd[v] - es[v].len <= budget
	 *
	 * We subtract es[v].len because v's sequence length is counted
	 * once in fwd[v] (as the path ...->v) and once in bwd[v]
	 * (as the path v->...), so it would be double-counted.
	 *
	 * We iterate the forward map and check each vertex against
	 * the backward map. kh_end(bwd) is the sentinel for "not found".
	 */
	for (k = 0; k < kh_end(fwd); ++k) {
		if (kh_exist(fwd, k)) {
			// kh_key works here because gwf_map64_t uses
			// KHASHL_MAP_INIT (buckets have .key/.val fields),
			// unlike gwf_set64_t which uses bare KHASHL_INIT.
			uint32_t v = (uint32_t)kh_key(fwd, k);
			int32_t df = kh_val(fwd, k);    // min dist v0->v
			khint_t kb = gwf_map64_get(bwd, (uint64_t)v);
			if (kb < kh_end(bwd)) {          // v reachable from v1
				int32_t db = kh_val(bwd, kb); // min dist v->v1
				if (df + db - es[v].len <= budget)
					gwf_set64_put(result, (uint64_t)v,
						&absent);
			}
		}
	}

	// Cleanup: free queue array and both hash maps
	kfree(km, queue.a);
	gwf_map64_destroy(fwd);
	gwf_map64_destroy(bwd);
}

typedef struct {
	const gfa_t *g; //graph
	const gfa_edseq_t *es; //edit sequence?
	const gfa_edopt_t *opt; //options
	int32_t ql; // query length
	const char *q; // query
	gwf_edbuf_t buf;
	int32_t s, n_a;
	gwf_diag_t *a; //oh, is that what a's are? diagonals?
	int32_t end_tb;
	uint32_t v0;
} gfa_edbuf_t;

void *gfa_ed_init(void *km, const gfa_edopt_t *opt, const gfa_t *g, const gfa_edseq_t *es, int32_t ql, const char *q, uint32_t v0, int32_t off0)
{
	gfa_edbuf_t *z;
	KCALLOC(km, z, 1);
	z->buf.km = km;
	z->opt = opt;
	z->g = g, z->es = es;
	z->ql = ql, z->q = q;
	z->buf.ha = gwf_set64_init2(km);
	z->buf.ht = gwf_map64_init2(km);
	kv_resize(gwf_trace_t, km, z->buf.t, 16);
	KCALLOC(km, z->a, 1);
	z->v0 = v0;
	z->a[0].vd = gwf_gen_vd(v0, -off0), z->a[0].k = off0 - 1, z->a[0].xo = 0;
	if (z->opt->traceback) z->a[0].t = gwf_trace_push(km, &z->buf.t, -1, -1, z->buf.ht);
	z->n_a = 1;
	return z;
}

/*
 * s_term seems to be the maximum edit distance allowed in the dynamic
 *        progrmming matrix. It is a termination condition which triggers if the
 *        score (s) in H_{s,v,k} becomes too high
 * i_term this is within z->opt. It looks like another early stopping condition
 */
void gfa_ed_step(void *z_, uint32_t v1, int32_t off1, int32_t s_term, gfa_edrst_t *r)
{
	FILE *stats_fp = gwf_ed_stats_fp();
	gfa_edbuf_t *z = (gfa_edbuf_t*)z_;
	const gfa_edopt_t *opt = z->opt;
	if (s_term < 0 && z->opt->s_term >= 0) s_term = z->opt->s_term;
	s_term = 3000;
	r->n_end = 0, r->n_iter = 0;

	// BFS reachable subgraph — computed before the wavefront
	// loop so we can assert the wavefront is a subset.
	gwf_set64_t *bfs_set = NULL;
	clock_t clk0 = 0, clk1 = 0;
	if (stats_fp && v1 != (uint32_t)-1) {
		clk0 = clock();
		bfs_set = gwf_set64_init2(z->buf.km);
		gwf_bfs_reachable(z->buf.km, z->g, z->es,
			z->v0, v1, z->ql + s_term, bfs_set);
		clk1 = clock();
	}

  //This for loop consists of the extend funciton, and a bunch of termination
  //condition behaviour
	while (z->n_a > 0) {
		z->a = gwf_ed_extend(&z->buf, opt, z->g, z->es, z->s, z->ql, z->q, v1, off1, &z->end_tb, &z->n_a, z->a, r);
		r->n_iter += z->n_a; // + z->buf.intv.n;
		if (r->end_off >= 0 || z->n_a == 0) break; //end offset probably only 
                                               //becomes >=0 when you reach the
                                               //end of a node. So this is
                                               //standard stopping condition
		if (r->n_end > 0) break;
		if (s_term >= 0 && z->s >= s_term) break; //terminate if score is outta ctrl
		if (z->opt->i_term > 0 && r->n_iter > z->opt->i_term) break;
		++z->s; //I guess this is (z->s)++. It makes me think s is the score
		if (gfa_ed_dbg >= 1) {
			printf("[%s] dist=%d, n=%d, n_intv=%ld, n_tb=%ld\n", __func__, z->s, z->n_a, z->buf.intv.n, z->buf.t.n);
			if (gfa_ed_dbg == 2) gwf_ed_print_diag(z->g, z->n_a, z->a);
			if (gfa_ed_dbg == 3) gwf_ed_print_intv(z->buf.intv.n, z->buf.intv.a);
		}
	}
	if (opt->traceback && r->end_off >= 0)
		gwf_traceback(&z->buf, r->end_v, z->end_tb, r);
	r->s = r->end_v != (uint32_t)-1? z->s : -1;

	// Print subgraph statistics
	if (stats_fp) {
		gwf_set64_t *hv, *hs;
		khint_t k;
		int absent;
		int32_t n_vtx = 0, n_seg = 0, n_edges = 0;
		int64_t bytes_full = 0, bytes_ess = 0, bytes_comp = 0;
		int64_t idx_arc_bytes = 0;

		// Extract unique vertices from buf->ha
		hv = gwf_set64_init2(z->buf.km);
		for (k = 0; k < kh_end(z->buf.ha); ++k) {
			if (kh_exist(z->buf.ha, k)) {
				uint32_t v = z->buf.ha->keys[k] >> 32;
				gwf_set64_put(hv, (uint64_t)v, &absent);
			}
		}

		// Assert: every wavefront vertex is BFS-reachable
		if (bfs_set) {
			for (k = 0; k < kh_end(hv); ++k) {
				if (kh_exist(hv, k)) {
					uint32_t v = (uint32_t)hv->keys[k];
					khint_t bk = gwf_set64_get(
						bfs_set, (uint64_t)v);
					assert(bk < kh_end(bfs_set));
				}
			}
		}

		// Build unique segment set from vertices
		hs = gwf_set64_init2(z->buf.km);
		for (k = 0; k < kh_end(hv); ++k) {
			if (kh_exist(hv, k)) {
				uint32_t v = (uint32_t)hv->keys[k];
				gwf_set64_put(hs, (uint64_t)(v >> 1), &absent);
			}
		}

		// Per-vertex costs: index + arcs
		for (k = 0; k < kh_end(hv); ++k) {
			if (kh_exist(hv, k)) {
				uint32_t v = (uint32_t)hv->keys[k];
				int32_t nv = gfa_arc_n(z->g, v);
				n_vtx++;
				n_edges += nv;
				idx_arc_bytes += sizeof(uint64_t)
					+ (int64_t)nv * sizeof(gfa_arc_t);
			}
		}

		// Per-segment costs: struct + sequence + name
		for (k = 0; k < kh_end(hs); ++k) {
			if (kh_exist(hs, k)) {
				uint32_t s = (uint32_t)hs->keys[k];
				int32_t slen = z->g->seg[s].len;
				n_seg++;
				bytes_full += sizeof(gfa_seg_t) + slen
					+ strlen(z->g->seg[s].name) + 1;
				bytes_ess += sizeof(int32_t) + slen;
				bytes_comp += sizeof(int32_t) + (slen + 3) / 4;
			}
		}
		bytes_full += idx_arc_bytes;
		bytes_ess += idx_arc_bytes;
		bytes_comp += idx_arc_bytes;

		fprintf(stats_fp, "TOTAL_SEGMENTS: %d\n", n_seg);
		fprintf(stats_fp, "TOTAL_VERTICES: %d\n", n_vtx);
		fprintf(stats_fp, "TOTAL_EDGES: %d\n", n_edges);
		fprintf(stats_fp, "BYTES_FULL: %lld\n",
			(long long)bytes_full);
		fprintf(stats_fp, "BYTES_ESSENTIAL: %lld\n",
			(long long)bytes_ess);
		fprintf(stats_fp, "BYTES_COMPRESSED: %lld\n",
			(long long)bytes_comp);
		fflush(stats_fp);

		gwf_set64_destroy(hs);
		gwf_set64_destroy(hv);

		/*
		 * BFS subgraph stats.
		 *
		 * bfs_set was built before the while loop. Here we
		 * just compute and print stats from it.
		 */
		if (bfs_set) {
			double bfs_ms;
			gwf_set64_t *bfs_seg;
			int32_t bfs_nv = 0, bfs_ns = 0, bfs_ne = 0;
			int64_t bfs_comp = 0, bfs_idx_arc = 0;
			int bfs_absent;

			bfs_ms = (double)(clk1 - clk0)
				/ CLOCKS_PER_SEC * 1000.0;

			/*
			 * Vertices are oriented (vertex_id = seg_id<<1|ori),
			 * so we deduplicate into segments (v>>1) to avoid
			 * counting both orientations of the same sequence.
			 *
			 * Compressed size =
			 *   per segment: sizeof(int32_t) + ceil(len/4)
			 *   per vertex: sizeof(uint64_t) for idx entry
			 *     + n_arcs * sizeof(gfa_arc_t)
			 */
			bfs_seg = gwf_set64_init2(z->buf.km);
			for (k = 0; k < kh_end(bfs_set); ++k) {
				if (kh_exist(bfs_set, k)) {
					uint32_t v = (uint32_t)bfs_set->keys[k];
					int32_t nv = gfa_arc_n(z->g, v);
					bfs_nv++;
					bfs_ne += nv;
					bfs_idx_arc += sizeof(uint64_t)
						+ (int64_t)nv * sizeof(gfa_arc_t);
					gwf_set64_put(bfs_seg,
						(uint64_t)(v >> 1), &bfs_absent);
				}
			}
			for (k = 0; k < kh_end(bfs_seg); ++k) {
				if (kh_exist(bfs_seg, k)) {
					uint32_t s = (uint32_t)bfs_seg->keys[k];
					int32_t slen = z->g->seg[s].len;
					bfs_ns++;
					bfs_comp += sizeof(int32_t)
						+ (slen + 3) / 4;
				}
			}
			bfs_comp += bfs_idx_arc;

			fprintf(stats_fp, "BFS_TIME_MS: %.3f\n", bfs_ms);
			fprintf(stats_fp, "BFS_VERTICES: %d\n", bfs_nv);
			fprintf(stats_fp, "BFS_SEGMENTS: %d\n", bfs_ns);
			fprintf(stats_fp, "BFS_EDGES: %d\n", bfs_ne);
			fprintf(stats_fp, "BFS_BYTES_COMPRESSED: %lld\n",
				(long long)bfs_comp);
			fflush(stats_fp);

			gwf_set64_destroy(bfs_seg);
		}
		if (bfs_set)
			gwf_set64_destroy(bfs_set);
	}
}

void gfa_ed_destroy(void *z_)
{
	gfa_edbuf_t *z = (gfa_edbuf_t*)z_;
	void *km = z->buf.km;
	kfree(km, z->a);
	gwf_set64_destroy(z->buf.ha);
	gwf_map64_destroy(z->buf.ht);
	kfree(km, z->buf.ooo.a);
	kfree(km, z->buf.intv.a);
	kfree(km, z->buf.tmp.a);
	kfree(km, z->buf.swap.a);
	kfree(km, z->buf.t.a);
	kfree(km, z);
}

int32_t gfa_edit_dist(void *km, const gfa_edopt_t *opt, const gfa_t *g, const gfa_edseq_t *es, int32_t ql, const char *q, uint32_t v0, int32_t off0, gfa_edrst_t *rst)
{
	void *z;
	z = gfa_ed_init(km, opt, g, es, ql, q, v0, off0);
	gfa_ed_step(z, (uint32_t)-1, -1, -1, rst);
	gfa_ed_destroy(z);
	return rst->s;
}
