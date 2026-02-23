#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>
#include "gfa-priv.h"
#include "kalloc.h"
#include "ksort.h"
#include "khashl.h" // make it compatible with kalloc
#include "kvec-km.h"

#ifndef GFA_ED_DBG
#define GFA_ED_DBG 0
#endif
int gfa_ed_dbg = GFA_ED_DBG;

/*--- gwfa input dump (opt-in via -DDUMP_GWFA) ---*/
#ifdef DUMP_GWFA
#include <sys/stat.h>

enum {
	D_QL, D_Q, D_STARTV, D_STARTOFF,
	D_ENDV, D_ENDOFF, D_STERM, D_DBG,
	D_NVTX, D_NARC, D_GRAPHSEQ,
	D_SEQOFF, D_SEQLEN,
	D_ARCV, D_ARCW, D_ARCOW, D_IDX,
	D_NFILES
};
static const char *dump_names[D_NFILES] = {
	"ql.txt", "q.txt", "startV.txt", "startOff.txt",
	"endV.txt", "endOff.txt", "s_term.txt", "dbg.txt",
	"n_vtx.txt", "n_arc.txt", "graphSeq.txt",
	"seq_off.txt", "seq_len.txt",
	"arc_v.txt", "arc_w.txt", "arc_ow.txt", "idx.txt"
};
static FILE *dump_fps[D_NFILES];
static int dump_inited;

static void dump_init(void) {
	int i;
	char path[512];
	if (dump_inited) return;
	mkdir("GwfaDump", 0755);
	for (i = 0; i < D_NFILES; i++) {
		snprintf(path, sizeof(path),
			"GwfaDump/%s", dump_names[i]);
		dump_fps[i] = fopen(path, "w");
	}
	dump_inited = 1;
}

static void dump_gwfa_inputs(
	int32_t ql, const char *q,
	uint32_t startV, int32_t startOff,
	uint32_t endV, int32_t endOff,
	subgfa_subgraph_t *sub,
	int32_t s_term, int dbg)
{
	uint32_t i;
	dump_init();
	/* scalars */
	fprintf(dump_fps[D_QL], "%d\n", ql);
	fprintf(dump_fps[D_Q], "%.*s\n", ql, q);
	fprintf(dump_fps[D_STARTV], "%u\n", startV);
	fprintf(dump_fps[D_STARTOFF], "%d\n", startOff);
	fprintf(dump_fps[D_ENDV], "%u\n", endV);
	fprintf(dump_fps[D_ENDOFF], "%d\n", endOff);
	fprintf(dump_fps[D_STERM], "%d\n", s_term);
	fprintf(dump_fps[D_DBG], "%d\n", dbg);
	/* subgraph */
	if (!sub) {
		fprintf(dump_fps[D_NVTX], "0\n");
		fprintf(dump_fps[D_NARC], "0\n");
		fprintf(dump_fps[D_GRAPHSEQ], "\n");
		fprintf(dump_fps[D_SEQOFF], "\n");
		fprintf(dump_fps[D_SEQLEN], "\n");
		fprintf(dump_fps[D_ARCV], "\n");
		fprintf(dump_fps[D_ARCW], "\n");
		fprintf(dump_fps[D_ARCOW], "\n");
		fprintf(dump_fps[D_IDX], "\n");
		return;
	}
	fprintf(dump_fps[D_NVTX], "%u\n", sub->n_vtx);
	fprintf(dump_fps[D_NARC],
		"%" PRIu64 "\n", sub->n_arc);
	/* graphSeq: compute total length from last vtx */
	{
		uint32_t total = 0;
		if (sub->n_vtx > 0)
			total = sub->seq_off[sub->n_vtx - 1]
				+ sub->seq_len[sub->n_vtx - 1];
		fprintf(dump_fps[D_GRAPHSEQ],
			"%.*s\n", (int)total, sub->graphSeq);
	}
	/* seq_off */
	for (i = 0; i < sub->n_vtx; i++)
		fprintf(dump_fps[D_SEQOFF], "%s%u",
			i ? " " : "", sub->seq_off[i]);
	fprintf(dump_fps[D_SEQOFF], "\n");
	/* seq_len */
	for (i = 0; i < sub->n_vtx; i++)
		fprintf(dump_fps[D_SEQLEN], "%s%d",
			i ? " " : "", sub->seq_len[i]);
	fprintf(dump_fps[D_SEQLEN], "\n");
	/* arcs */
	for (i = 0; i < (uint32_t)sub->n_arc; i++)
		fprintf(dump_fps[D_ARCV], "%s%u",
			i ? " " : "", sub->arc[i].v);
	fprintf(dump_fps[D_ARCV], "\n");
	for (i = 0; i < (uint32_t)sub->n_arc; i++)
		fprintf(dump_fps[D_ARCW], "%s%u",
			i ? " " : "", sub->arc[i].w);
	fprintf(dump_fps[D_ARCW], "\n");
	for (i = 0; i < (uint32_t)sub->n_arc; i++)
		fprintf(dump_fps[D_ARCOW], "%s%d",
			i ? " " : "", sub->arc[i].ow);
	fprintf(dump_fps[D_ARCOW], "\n");
	/* idx */
	for (i = 0; i < sub->n_vtx; i++)
		fprintf(dump_fps[D_IDX], "%s%" PRIu64,
			i ? " " : "", sub->idx[i]);
	fprintf(dump_fps[D_IDX], "\n");
}

void dump_gwfa_flush(void) {
	int i;
	for (i = 0; i < D_NFILES; i++)
		if (dump_fps[i]) {
			fclose(dump_fps[i]);
			dump_fps[i] = NULL;
		}
}
#endif /* DUMP_GWFA */

static FILE *gwf_scores_fp(void) {
	static FILE *fp;
	if (!fp) fp = fopen("scores.txt", "a");
	return fp;
}

static FILE *gwf_wf_debug_fp(void) {
	static FILE *fp;
	if (!fp) fp = fopen("wfDebug.txt", "w");
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
			t[s->len - j - 1] =
				gfa_comp_table[(uint8_t)s->seq[j]];
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

/* Static array capacity constants */
#define DIAG_CAP  (16 << 20)   /* 16M = 2^24 */
#define INTV_CAP  (1 << 21)    /* 2M (~32 MB) */
#define HA_BITS   22           /* 4M slots */
#define HA_CAP    (1 << HA_BITS)
#define HA_MASK   (HA_CAP - 1)
#define A_MASK    (DIAG_CAP - 1)

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

#define subgfa_arc_key(a) ((uint64_t)(a).v << 32)
KRADIX_SORT_INIT(subgfa_arc, subgfa_arc_t, subgfa_arc_key, 8)

#define subgfa_arc_n(s, v) \
	((uint32_t)(s)->idx[(v)])
#define subgfa_arc_a(s, v) \
	(&(s)->arc[(s)->idx[(v)]>>32])

static int gwf_intv_is_sorted(int32_t n_a,
	const gwf_intv_t *a)
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
static size_t gwf_intv_merge2(gwf_intv_t *a,
	size_t n_b, const gwf_intv_t *b,
	size_t n_c, const gwf_intv_t *c)
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
 * Observe that there are no arrays. This implies that
 * it must be a single diagonal not a wavefront
 */
typedef struct { // a diagonal
	uint64_t vd; // higher 32 bits: vertex ID; lower 32 bits: diagonal+0x4000000
	int32_t k; // wavefront position on the vertex
} gwf_diag_t;

typedef kvec_t(gwf_diag_t) gwf_diag_v;

#define ed_key(x) ((x).vd)
KRADIX_SORT_INIT(gwf_ed, gwf_diag_t, ed_key, 8)

KHASHL_MAP_INIT(KH_LOCAL, gwf_map64_t, gwf_map64,
	uint64_t, int32_t, kh_hash_uint64, kh_eq_generic)

/*
 * Core GWFA routine
 */
KHASHL_INIT(KH_LOCAL, gwf_set64_t, gwf_set64,
	uint64_t, kh_hash_dummy, kh_eq_generic)

/* File-scope static arrays */
static gwf_diag_t s_diag_a[DIAG_CAP];
static gwf_diag_t s_diag_b[DIAG_CAP];
static gwf_diag_t s_A[DIAG_CAP];
static gwf_intv_t s_intv[INTV_CAP];
static gwf_intv_t s_tmp[INTV_CAP];
static gwf_intv_t s_swap[INTV_CAP];
static gwf_diag_t s_sort_buf[DIAG_CAP];
static uint64_t   s_ha_keys[HA_CAP];
static uint8_t    s_ha_occ[HA_CAP];
static uint32_t   s_ha_dirty[HA_CAP];

/* File-scope mutable counters */
static uint32_t s_ha_n_dirty;
static uint32_t s_A_head, s_A_tail, s_A_count;
static size_t   s_intv_n, s_tmp_n;

/* ---- hash set helpers ---- */
static inline void ha_clear(void) {
	uint32_t i;
	for (i = 0; i < s_ha_n_dirty; ++i)
		s_ha_occ[s_ha_dirty[i]] = 0;
	s_ha_n_dirty = 0;
}

static inline uint32_t ha_put(
	uint64_t key, int *absent)
{
	uint32_t h = (uint32_t)key * 2654435769U
		>> (32 - HA_BITS);
	while (s_ha_occ[h]) {
		if (s_ha_keys[h] == key) {
			*absent = 0;
			return h;
		}
		h = (h + 1) & HA_MASK;
	}
	s_ha_keys[h] = key;
	s_ha_occ[h] = 1;
	s_ha_dirty[s_ha_n_dirty++] = h;
	*absent = 1;
	return h;
}

/* ---- queue A helpers ---- */
static inline void A_clear(void) {
	s_A_head = s_A_tail = s_A_count = 0;
}

static inline uint32_t A_size(void) {
	return s_A_count;
}

static inline gwf_diag_t *A_pushp(void) {
	gwf_diag_t *p =
		&s_A[s_A_tail++ & A_MASK];
	s_A_count++;
	return p;
}

static inline gwf_diag_t A_shift(void) {
	gwf_diag_t t =
		s_A[s_A_head++ & A_MASK];
	s_A_count--;
	return t;
}

// push (v,d,k) to the end of an array
static inline void gwf_diag_push(
	gwf_diag_v *a, uint32_t v, int32_t d, int32_t k)
{
	gwf_diag_t *p = &a->a[a->n++];
	p->vd = gwf_gen_vd(v, d), p->k = k;
}

// determine the wavefront on diagonal (v,d)
static inline int32_t gwf_diag_update(gwf_diag_t *p,
	uint32_t v, int32_t d, int32_t k)
{
	uint64_t vd = gwf_gen_vd(v, d);
	if (p->vd == vd) {
		p->k = p->k > k ? p->k : k;
		return 0;
	}
	return 1;
}

// sort a[] using n_sorted as the sorted/unsorted split point.
static void gwf_diag_sort(int32_t n_a, gwf_diag_t *a,
	int32_t n_sorted, gwf_diag_v *buf)
{
	int32_t i, j, k, n_b, n_c;
	gwf_diag_t *b, *c;

	n_b = n_sorted;
	n_c = n_a - n_sorted;
	b = buf->a, c = b + n_b;
	memcpy(b, a, n_b * sizeof(*a));
	memcpy(c, a + n_b, n_c * sizeof(*a));
	radix_sort_gwf_ed(c, c + n_c);

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
static int32_t gwf_diag_dedup(int32_t n_a,
	gwf_diag_t *a, int32_t n_sorted,
	gwf_diag_v *buf)
{
	int32_t i, n, st;
	if (n_sorted < n_a)
		gwf_diag_sort(n_a, a, n_sorted, buf);
	for (i = 1, st = 0, n = 0; i <= n_a; ++i) {
		if (i == n_a || a[i].vd != a[st].vd) {
			int32_t j, max_j = st;
			if (st + 1 < i)
				for (j = st + 1; j < i; ++j) // choose the far end (i.e. the wavefront)
					if (a[max_j].k < a[j].k)
						max_j = j;
			a[n++] = a[max_j];
			st = i;
		}
	}
	return n;
}

// use forbidden bands to remove diagonals not on the wavefront
static int32_t gwf_mixed_dedup(int32_t n_a,
	gwf_diag_t *a, int32_t n_b, gwf_intv_t *b)
{
	int32_t i = 0, j = 0, k = 0;
	while (i < n_a && j < n_b) {
		if (a[i].vd >= b[j].vd0
			&& a[i].vd < b[j].vd1) ++i;
		else if (a[i].vd >= b[j].vd1) ++j;
		else a[k++] = a[i++];
	}
	while (i < n_a) a[k++] = a[i++];
	return k;
}


// remove diagonals not on the wavefront
static int32_t gwf_dedup(int32_t n_a,
	gwf_diag_t *a, int32_t n_sorted)
{
	if (s_intv_n + s_tmp_n > 0) {
		size_t swap_n;
		if (!gwf_intv_is_sorted(s_tmp_n, s_tmp))
			radix_sort_gwf_intv(s_tmp,
				s_tmp + s_tmp_n);
		memcpy(s_swap, s_intv,
			s_intv_n * sizeof(gwf_intv_t));
		swap_n = s_intv_n;
		s_intv_n = gwf_intv_merge2(s_intv,
			swap_n, s_swap, s_tmp_n, s_tmp);
	}
	gwf_diag_v sb;
	sb.a = s_sort_buf; sb.n = 0; sb.m = DIAG_CAP;
	n_a = gwf_diag_dedup(n_a, a, n_sorted, &sb);
	if (s_intv_n > 0)
		n_a = gwf_mixed_dedup(n_a, a,
			s_intv_n, s_intv);
	return n_a;
}


//zkn this is the function that pushes the wavefront
//along the diagonal with matches.
//We do examine the query here
// reach the wavefront
static inline int32_t gwf_extend1(int32_t d, int32_t k,
	int32_t vl, const char *ts,
	int32_t ql, const char *qs)
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
static void gwf_ed_extend_batch(
	const subgfa_subgraph_t *sub,
	int32_t ql, const char *q, int32_t n,
	gwf_diag_t *a, gwf_diag_v *B)
{
	int32_t j, m;
	int32_t v = a->vd>>32;
	int32_t vl = sub->seq_len[v];
	const char *ts = sub->graphSeq + sub->seq_off[v];
	gwf_diag_t *b;

	// wfa_extend
	for (j = 0; j < n; ++j)
		a[j].k = gwf_extend1((int32_t)a[j].vd
			- GWF_DIAG_SHIFT, a[j].k, vl, ts,
			ql, q);

	// wfa_next
	b = &B->a[B->n];
	b[0].vd = a[0].vd - 1;
	b[0].k = a[0].k + 1;
	b[1].vd = a[0].vd;
	b[1].k = (n == 1 || a[0].k > a[1].k
		? a[0].k : a[1].k) + 1;
	for (j = 1; j < n - 1; ++j) {
		int32_t k = a[j-1].k;
		k = k > a[j].k + 1 ? k : a[j].k + 1;
		k = k > a[j+1].k + 1 ? k : a[j+1].k + 1;
		b[j+1].vd = a[j].vd, b[j+1].k = k;
	}
	if (n >= 2) {
		b[n].vd = a[n-1].vd;
		b[n].k = a[n-2].k > a[n-1].k + 1
			? a[n-2].k : a[n-1].k + 1;
	}
	b[n+1].vd = a[n-1].vd + 1;
	b[n+1].k = a[n-1].k;

	// drop out-of-bound cells
	//if (a[n-1].k == vl - 1) b[n+1].k = vl; // insertion to the end of a vertex is handled elsewhere. FIXME: this line leads to wrong result for MHC-57 and MHC-HG002.2
	for (j = 0; j < n; ++j) {
		gwf_diag_t *p = &a[j];
		if (p->k == vl - 1 || (int32_t)p->vd
			- GWF_DIAG_SHIFT + p->k == ql - 1)
			*A_pushp() = *p;
	}
	for (j = 0, m = 0; j < n + 2; ++j) {
		gwf_diag_t *p = &b[j];
		int32_t d = (int32_t)p->vd - GWF_DIAG_SHIFT;
		if (d + p->k < ql && p->k < vl) {
			b[m++] = *p;
		} else if (p->k == vl) {
			gwf_intv_t *qi;
			qi = &s_tmp[s_tmp_n++];
			qi->vd0 = gwf_gen_vd(v, d);
			qi->vd1 = qi->vd0 + 1;
		}
	}
	B->n += m;
}

//This function demonstrates, 1) that we are not doing
//affine gap, and 2) that we are doing base level
//alignment
// wfa_extend and wfa_next combined
static gwf_diag_t *gwf_ed_extend(
	const subgfa_subgraph_t *sub,
	int32_t s, int32_t ql, const char *q,
	uint32_t endV, int32_t endOff, int32_t *n_a_,
	gwf_diag_t *a, int *terminate)
{
	int32_t i, x, n = *n_a_, do_dedup = 1;
	gwf_diag_v B;
	gwf_diag_t *b;

	s_tmp_n = 0;
	ha_clear(); // hash table $h to avoid visiting a vertex twice
	A_clear(); // $A is a queue

	/* B is the OTHER static buffer (ping-pong) */
	B.a = (a == s_diag_a) ? s_diag_b : s_diag_a;
	B.n = 0;
	B.m = DIAG_CAP;

	for (x = 0, i = 1; i <= n; ++i) {
		if (i == n || a[i].vd != a[i-1].vd + 1) {
			gwf_ed_extend_batch(sub, ql, q,
				i - x, &a[x], &B);
			x = i;
		}
	}
	if (A_size() == 0) do_dedup = 0;
	int32_t n_sorted = B.n;
	// $a is not used as it has been copied to $A

	while (A_size()) {
		gwf_diag_t t;
		uint32_t v;
		int32_t d, k, i, vl;

		t = A_shift();
		v = t.vd >> 32; // vertex
		d = (int32_t)t.vd - GWF_DIAG_SHIFT; // diagonal
		k = t.k; // wavefront position on the vertex
		vl = sub->seq_len[v];
		k = gwf_extend1(d, k, vl,
			sub->graphSeq + sub->seq_off[v],
			ql, q);
		i = k + d; // query position

		if (k + 1 < vl && i + 1 < ql) { // the most common case: the wavefront is in the middle
			int32_t push1 = 1, push2 = 1;
			if (B.n >= 2) push1 = gwf_diag_update(
				&B.a[B.n - 2], v, d-1, k+1);
			if (B.n >= 1) push2 = gwf_diag_update(
				&B.a[B.n - 1], v, d,   k+1);
			if (push1)
				gwf_diag_push(&B, v, d-1, k+1);
			if (push2 || push1)
				gwf_diag_push(&B, v, d,   k+1);
			gwf_diag_push(&B, v, d+1, k);
    //This else is the branchy weirdness part where you keep pushing along edges
		} else if (i + 1 < ql) { // k + 1 == vl; reaching the end of the vertex but not the end of query
			int32_t nv = subgfa_arc_n(sub, v);
			int32_t j, n_ext = 0;
			const subgfa_arc_t *av =
				subgfa_arc_a(sub, v);
			gwf_intv_t *p;
			p = &s_tmp[s_tmp_n++];
			p->vd0 = gwf_gen_vd(v, d);
			p->vd1 = p->vd0 + 1;
			for (j = 0; j < nv; ++j) {
				uint32_t w = av[j].w;
				int32_t ol = av[j].ow;
				int absent;
				ha_put(
					(uint64_t)w<<32 | (i + 1),
					&absent);
				if (q[i + 1] == (sub->graphSeq
					+ sub->seq_off[w])[ol]) {
					++n_ext;
					if (absent) {
						gwf_diag_t *dp;
						dp = A_pushp();
						dp->vd = gwf_gen_vd(
							w, i + 1 - ol);
						dp->k = ol;
					}
				} else if (absent) {
					gwf_diag_push(&B,
						w, i - ol, ol);
					gwf_diag_push(&B,
						w, i + 1 - ol, ol);
				}
			}
			if (nv == 0 || n_ext != nv)
				gwf_diag_push(&B,
					v, d+1, k);
		} else if (endV == (uint32_t)-1
			|| (v == endV && k == endOff)) { // i + 1 == ql
			*terminate = 1;
			return 0;
		} else if (k + 1 < vl) { // i + 1 == ql; reaching the end of the query but not the end of the vertex
			gwf_diag_push(&B, v, d-1, k+1);
		} else if (v != endV) { // i + 1 == ql && k + 1 == vl; not reaching the last vertex $endV
			int32_t nv = subgfa_arc_n(sub, v), j;
			const subgfa_arc_t *av =
				subgfa_arc_a(sub, v);
			for (j = 0; j < nv; ++j)
				gwf_diag_push(&B,
					av[j].w, i - av[j].ow,
					av[j].ow);
		} else { // may come here when k>endOff (due to banding); do nothing in this case
		}
	}

	*n_a_ = n = B.n, b = B.a;

	if (do_dedup)
		*n_a_ = n = gwf_dedup(n, b, n_sorted);
	return b;
}


static void gwf_ed_print_diag(const gfa_t *g,
	size_t n, gwf_diag_t *a) // for debugging only
{
	size_t i;
	for (i = 0; i < n; ++i) {
		int32_t d = (int32_t)a[i].vd
			- GWF_DIAG_SHIFT;
		printf("Z\t%d\t%s\t%d\t%d\n", d + a[i].k,
			g->seg[(a[i].vd>>32)>>1].name,
			d, a[i].k);
	}
}

static void gwf_ed_print_intv(size_t n,
	gwf_intv_t *a) // for debugging only
{
	FILE *fp = gwf_wf_debug_fp();
	size_t i;
	for (i = 0; i < n; ++i)
		fprintf(fp, "Z\t%d\t%d\t%d\n",
			(int32_t)(a[i].vd0>>32),
			(int32_t)a[i].vd0 - GWF_DIAG_SHIFT,
			(int32_t)a[i].vd1 - GWF_DIAG_SHIFT);
	fflush(fp);
}

static void gwf_ed_print_wf(int32_t n,
	const gwf_diag_t *a)
{
	FILE *fp = gwf_wf_debug_fp();
	int32_t i;
	for (i = 0; i < n; ++i) {
		int32_t nid = (int32_t)(a[i].vd >> 32);
		int32_t diag = (int32_t)a[i].vd
			- GWF_DIAG_SHIFT;
		fprintf(fp, "WF\t%d\t%d\t%d\n",
			nid, diag, a[i].k);
	}
	fflush(fp);
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
 *     NOTE* BFS, but we're doing unweighted. so we still
 *     need to check when we revist a node it could be
 *     shorter
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
 * @param km      thread-local memory allocator (NULL ok)
 * @param g       the graph (for arc traversal)
 * @param es      per-vertex sequences (for lengths)
 * @param v0      source vertex (start of GWFA alignment)
 * @param v1      sink vertex (target of GWFA alignment)
 * @param budget  max total sequence length (= query length)
 * @param result  output: set of vertices on valid paths
 */
static void gwf_bfs_reachable(void *km,
	const gfa_t *g, const gfa_edseq_t *es,
	uint32_t v0, uint32_t v1,
	int32_t off0, int32_t off1,
	int32_t budget, gwf_set64_t *result)
{
	gwf_map64_t *fwd, *bwd;
	typedef kvec_t(uint32_t) uint32_v;
	uint32_v queue = {0, 0, 0};
	int32_t head, j;
	khint_t k;
	int absent;

	/*
	 * PASS 1: Forward BFS from v0.
	 *
	 * Seed v0 with its effective length (es[v0].len - off0)
	 * since the wavefront starts at offset off0 inside v0.
	 * Don't prune v1: its full-segment fwd distance can
	 * exceed budget, but the intersection subtracts the
	 * excess via bwd[v1] = off1 + 1 (< es[v1].len).
	 */
	fwd = gwf_map64_init2(km);
	k = gwf_map64_put(fwd, (uint64_t)v0, &absent);
	kh_val(fwd, k) = es[v0].len - off0;
	kv_push(uint32_t, km, queue, v0);

	for (head = 0; head < (int32_t)queue.n; ++head) {
		uint32_t v = queue.a[head];
		khint_t kv = gwf_map64_get(fwd,
			(uint64_t)v);
		int32_t dv = kh_val(fwd, kv);
		int32_t nv = gfa_arc_n(g, v);
		const gfa_arc_t *av = gfa_arc_a(g, v);
		for (j = 0; j < nv; ++j) {
			uint32_t w = av[j].w;
			int32_t dw = dv + es[w].len;
			khint_t kw;
			// Don't prune v1: it may exceed budget in
			// fwd but the intersection corrects for it.
			if (dw > budget && w != v1) continue;
			kw = gwf_map64_put(fwd, (uint64_t)w,
				&absent);
			if (absent || kh_val(fwd, kw) > dw) {
				kh_val(fwd, kw) = dw;
				kv_push(uint32_t, km, queue, w);
			}
		}
	}

	/*
	 * PASS 2: Backward BFS from v1 (reverse edges).
	 *
	 * Seed v1 with off1 + 1 since the wavefront ends at
	 * offset off1 inside v1. Don't prune v0: same reason
	 * as above (the intersection corrects for it).
	 */
	bwd = gwf_map64_init2(km);
	k = gwf_map64_put(bwd, (uint64_t)v1, &absent);
	kh_val(bwd, k) = off1 + 1;
	queue.n = 0;
	kv_push(uint32_t, km, queue, v1);

	for (head = 0; head < (int32_t)queue.n; ++head) {
		uint32_t v = queue.a[head];
		khint_t kv = gwf_map64_get(bwd,
			(uint64_t)v);
		int32_t dv = kh_val(bwd, kv);
		int32_t nv = gfa_arc_n(g, v ^ 1);
		const gfa_arc_t *av = gfa_arc_a(g, v ^ 1);
		for (j = 0; j < nv; ++j) {
			uint32_t pred = av[j].w ^ 1;
			int32_t dp = dv + es[pred].len;
			khint_t kp;
			// Don't prune v0: same logic as v1 above.
			if (dp > budget && pred != v0)
				continue;
			kp = gwf_map64_put(bwd,
				(uint64_t)pred, &absent);
			if (absent || kh_val(bwd, kp) > dp) {
				kh_val(bwd, kp) = dp;
				kv_push(uint32_t, km, queue,
					pred);
			}
		}
	}

	/*
	 * INTERSECTION: vertex v is on a valid v0->v1 path iff
	 *   fwd[v] + bwd[v] - es[v].len <= budget
	 *
	 * With adjusted seeds this gives the true path length
	 * from off0 in v0 to off1 in v1 through v.
	 */
	for (k = 0; k < kh_end(fwd); ++k) {
		if (kh_exist(fwd, k)) {
			uint32_t v = (uint32_t)kh_key(fwd, k);
			int32_t df = kh_val(fwd, k);
			khint_t kb = gwf_map64_get(bwd,
				(uint64_t)v);
			if (kb < kh_end(bwd)) {
				int32_t db = kh_val(bwd, kb);
				if (df + db - es[v].len <= budget)
					gwf_set64_put(result,
						(uint64_t)v, &absent);
			}
		}
	}

	{
		khint_t k0 = gwf_set64_get(result,
			(uint64_t)v0);
		khint_t k1 = gwf_set64_get(result,
			(uint64_t)v1);
		if (k0 >= kh_end(result)
			|| k1 >= kh_end(result)) {
			fprintf(stderr,
				"[gwf_bfs] v0=%u %s subgraph, "
				"v1=%u %s subgraph, budget=%d\n",
				v0,
				k0 < kh_end(result)
					? "IN" : "NOT IN",
				v1,
				k1 < kh_end(result)
					? "IN" : "NOT IN",
				budget);
			kfree(km, queue.a);
			gwf_map64_destroy(fwd);
			gwf_map64_destroy(bwd);
			return;
		}
	}

	kfree(km, queue.a);
	gwf_map64_destroy(fwd);
	gwf_map64_destroy(bwd);
}

/*
 * Build a compact subgraph from BFS-reachable vertices.
 */
subgfa_subgraph_t *subgfa_subgraph(const gfa_t *g,
	const gfa_edseq_t *es, uint32_t v0, uint32_t v1,
	int32_t off0, int32_t off1, int32_t budget,
	int32_t **seg_remap_out)
{
	void *km = NULL;
	gwf_set64_t *vset, *seg_set;
	subgfa_subgraph_t *sub;
	uint32_t *segs, n_seg, i;
	int32_t *seg_remap;
	uint64_t n_arc, arc_idx;
	khint_t k;
	int absent;

	// 1. BFS to get vertex hash set
	vset = gwf_set64_init2(km);
	gwf_bfs_reachable(km, g, es, v0, v1, off0, off1,
		budget, vset);

	if (kh_size(vset) == 0) {
		gwf_set64_destroy(vset);
		*seg_remap_out = NULL;
		return NULL;
	}

	// 2. Collect unique segment IDs
	segs = (uint32_t*)malloc(
		kh_size(vset) * sizeof(uint32_t));
	seg_set = gwf_set64_init2(km);
	n_seg = 0;
	for (k = 0; k < kh_end(vset); ++k) {
		if (kh_exist(vset, k)) {
			uint32_t seg = (uint32_t)vset->keys[k]
				>> 1;
			gwf_set64_put(seg_set, (uint64_t)seg,
				&absent);
			if (absent) segs[n_seg++] = seg;
		}
	}
	gwf_set64_destroy(seg_set);

	// 3. Sort segment IDs (insertion sort, n_seg small)
	for (i = 1; i < n_seg; ++i) {
		uint32_t tmp = segs[i], j = i;
		while (j > 0 && segs[j - 1] > tmp) {
			segs[j] = segs[j - 1];
			--j;
		}
		segs[j] = tmp;
	}

	// Build reverse mapping: old_seg -> new_seg
	seg_remap = (int32_t*)calloc(
		g->n_seg, sizeof(int32_t));
	for (i = 0; i < (uint32_t)g->n_seg; ++i)
		seg_remap[i] = -1;
	for (i = 0; i < n_seg; ++i)
		seg_remap[segs[i]] = (int32_t)i;

	// 4. Allocate subgraph
	sub = (subgfa_subgraph_t*)calloc(
		1, sizeof(subgfa_subgraph_t));
	sub->n_vtx = n_seg * 2;

	// 5. Build concatenated graphSeq (both strands)
	GFA_MALLOC(sub->seq_off, n_seg * 2);
	GFA_MALLOC(sub->seq_len, n_seg * 2);
	{
		uint32_t total_len = 0;
		for (i = 0; i < n_seg; ++i) {
			uint32_t old_s = segs[i];
			sub->seq_len[i<<1] =
				es[old_s<<1].len;
			sub->seq_len[i<<1|1] =
				es[old_s<<1|1].len;
			total_len += es[old_s<<1].len
				+ es[old_s<<1|1].len;
		}
		GFA_MALLOC(sub->graphSeq, total_len);
		total_len = 0;
		for (i = 0; i < n_seg; ++i) {
			uint32_t old_s = segs[i];
			int32_t len;
			len = es[old_s<<1].len;
			sub->seq_off[i<<1] = total_len;
			memcpy(sub->graphSeq + total_len,
				es[old_s<<1].seq, len);
			total_len += len;
			len = es[old_s<<1|1].len;
			sub->seq_off[i<<1|1] = total_len;
			memcpy(sub->graphSeq + total_len,
				es[old_s<<1|1].seq, len);
			total_len += len;
		}
	}

	// 7. Count arcs (both endpoints in vertex set)
	n_arc = 0;
	for (k = 0; k < kh_end(vset); ++k) {
		if (kh_exist(vset, k)) {
			uint32_t v = (uint32_t)vset->keys[k];
			int32_t nv = gfa_arc_n(g, v), j;
			const gfa_arc_t *av = gfa_arc_a(g, v);
			for (j = 0; j < nv; ++j) {
				khint_t kw = gwf_set64_get(
					vset, (uint64_t)av[j].w);
				if (kw < kh_end(vset))
					++n_arc;
			}
		}
	}
	sub->n_arc = n_arc;

	// 8. Copy arcs with remapped segment IDs
	sub->arc = (subgfa_arc_t*)malloc(
		n_arc * sizeof(subgfa_arc_t));
	arc_idx = 0;
	for (k = 0; k < kh_end(vset); ++k) {
		if (kh_exist(vset, k)) {
			uint32_t v = (uint32_t)vset->keys[k];
			int32_t nv = gfa_arc_n(g, v), j;
			const gfa_arc_t *av = gfa_arc_a(g, v);
			for (j = 0; j < nv; ++j) {
				khint_t kw = gwf_set64_get(
					vset, (uint64_t)av[j].w);
				if (kw < kh_end(vset)) {
					sub->arc[arc_idx].v =
						(seg_remap[v >> 1] << 1)
						| (v & 1);
					sub->arc[arc_idx].w =
						(seg_remap[av[j].w >> 1]
						<< 1)
						| (av[j].w & 1);
					sub->arc[arc_idx].ow =
						av[j].ow;
					++arc_idx;
				}
			}
		}
	}

	// 9. Sort arcs by source v
	radix_sort_subgfa_arc(sub->arc,
		sub->arc + n_arc);

	// 10. Build idx[] (size n_vtx)
	sub->idx = (uint64_t*)calloc(
		n_seg * 2, sizeof(uint64_t));
	if (n_arc > 0) {
		uint32_t cur_v = sub->arc[0].v;
		uint64_t start = 0;
		for (i = 1; i <= (uint32_t)n_arc; ++i) {
			if (i == (uint32_t)n_arc
				|| sub->arc[i].v != cur_v) {
				sub->idx[cur_v] =
					(start << 32)
					| (i - start);
				if (i < (uint32_t)n_arc) {
					cur_v = sub->arc[i].v;
					start = i;
				}
			}
		}
	}

	// 11. Free BFS temporaries; return seg_remap
	free(segs);
	*seg_remap_out = seg_remap;
	gwf_set64_destroy(vset);
	return sub;
}

void subgfa_subgraph_destroy(subgfa_subgraph_t *sub)
{
	if (!sub) return;
	free(sub->graphSeq);
	free(sub->seq_off);
	free(sub->seq_len);
	free(sub->arc);
	free(sub->idx);
	free(sub);
}

typedef struct {
	const gfa_t *g; //graph
	const gfa_edseq_t *es; //edit sequence?
	const gfa_edopt_t *opt; //options
	int32_t ql; // query length
	const char *q; // query
	uint32_t v0;
	int32_t off0;
} gfa_edbuf_t;

void *gfa_ed_init(void *km, const gfa_edopt_t *opt,
	const gfa_t *g, const gfa_edseq_t *es,
	int32_t ql, const char *q,
	uint32_t v0, int32_t off0)
{
	gfa_edbuf_t *z;
	(void)km;
	z = (gfa_edbuf_t*)calloc(1, sizeof(*z));
	z->opt = opt;
	z->g = g, z->es = es;
	z->ql = ql, z->q = q;
	z->v0 = v0;
	z->off0 = off0;
	return z;
}

int gwfa(int32_t ql, const char *q,
	uint32_t startV, int32_t startOff,
	uint32_t endV, int32_t endOff,
	subgfa_subgraph_t *sub, int32_t s_term,
	int dbg)
{
	gwf_diag_t *a;
	int32_t n_a, s;
	int terminate = 0;

	/* Reset per-alignment counters */
	s_intv_n = 0;
	s_tmp_n = 0;
	memset(s_ha_occ, 0, sizeof(s_ha_occ));
	s_ha_n_dirty = 0;

	/* Initial wavefront */
	a = s_diag_a;
	n_a = 1;
	a[0].vd = gwf_gen_vd(startV, -startOff);
	a[0].k = startOff - 1;

	s = 0;
	while (n_a > 0) {
		a = gwf_ed_extend(sub, s, ql, q,
			endV, endOff, &n_a, a, &terminate);
		if (terminate || s >= s_term) break;
		++s;
		if (dbg >= 1) {
			FILE *fp = gwf_wf_debug_fp();
			fprintf(fp,
				"[gfa_ed_step] dist=%d, n=%d,"
				" n_intv=%zd, n_tb=0\n",
				s, n_a, s_intv_n);
			fflush(fp);
			gwf_ed_print_intv(s_intv_n, s_intv);
			gwf_ed_print_wf(n_a, a);
		}
	}
	return terminate ? s : -1;
}

/*
 * s_term seems to be the maximum edit distance allowed in
 *        the dynamic programming matrix. It is a
 *        termination condition which triggers if the score
 *        (s) in H_{s,v,k} becomes too high
 * i_term this is within z->opt. It looks like another
 *        early stopping condition
 */
void gfa_ed_step(void *z_, uint32_t v1,
	int32_t off1, int32_t s_term, gfa_edrst_t *r)
{
	gfa_edbuf_t *z = (gfa_edbuf_t*)z_;
	int32_t *seg_remap = NULL;
	subgfa_subgraph_t *sub;
	uint32_t rv0, rv1;

	if (s_term < 0 && z->opt->s_term >= 0)
		s_term = z->opt->s_term;
	s_term = 3000;
	r->n_end = 0, r->n_iter = 0;

	assert(v1 != (uint32_t)-1);
	sub = subgfa_subgraph(z->g, z->es,
		z->v0, v1, z->off0, off1,
		z->ql + s_term, &seg_remap);

	if (sub == NULL) {
		r->s = -1;
		r->end_v = (uint32_t)-1;
		r->end_off = -1;
#ifdef DUMP_GWFA
		dump_gwfa_inputs(z->ql, z->q,
			0, z->off0, 0, off1,
			NULL, s_term, gfa_ed_dbg);
#endif
		FILE *fp = gwf_scores_fp();
		if (fp) {
			fprintf(fp, "%d\n", -1);
			fflush(fp);
		}
		return;
	}

	// Remap start/end vertices to compact vertex space
	rv0 = (seg_remap[z->v0 >> 1] << 1)
		| (z->v0 & 1);
	rv1 = (seg_remap[v1 >> 1] << 1)
		| (v1 & 1);

#ifdef DUMP_GWFA
	dump_gwfa_inputs(z->ql, z->q,
		rv0, z->off0, rv1, off1,
		sub, s_term, gfa_ed_dbg);
#endif
	int score = gwfa(z->ql, z->q,
		rv0, z->off0, rv1, off1,
		sub, s_term, gfa_ed_dbg);
	r->s = score;

	{
		FILE *fp = gwf_scores_fp();
		if (fp) {
			fprintf(fp, "%d\n", score);
			fflush(fp);
		}
	}

	subgfa_subgraph_destroy(sub);
	free(seg_remap);
}

void gfa_ed_destroy(void *z_)
{
	free(z_);
}

int32_t gfa_edit_dist(void *km,
	const gfa_edopt_t *opt, const gfa_t *g,
	const gfa_edseq_t *es, int32_t ql,
	const char *q, uint32_t v0, int32_t off0,
	gfa_edrst_t *rst)
{
	void *z;
	z = gfa_ed_init(km, opt, g, es, ql, q, v0, off0);
	gfa_ed_step(z, (uint32_t)-1, -1, -1, rst);
	gfa_ed_destroy(z);
	return rst->s;
}
