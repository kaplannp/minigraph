#ifndef SERIALIZATION_HH
#define SERIALIZATION

#include "minigraph.h"
#include "gfa-priv.h"

typedef struct {
	void *km; //this is a kalloc handler? 
	const gfa_t *g;
	const gfa_edseq_t *es; //A gfa_edseq_t is a sequence. literally char*+len.
                         //It is possible that this es is an array, I don't know
                         //what the count would be
	const char *qseq; //query sequence
	int32_t n_seg, n_llc, m_llc, n_a; //presumably lc is linear chain, so perahps 
                                    //this relates to the two chains we pass as
                                    //other arguments to the kernel
	mg_llchain_t *llc; //offset, count, vertex, score, edit distance. Note this is
                     //probably an array, therefore I think collectively these
                     //llc form a complete path through the graph
} bridge_aux_t;

#ifdef __cplusplus
#include <cstdint>
extern "C" {
#endif
void print( bridge_aux_t* aux, int32_t kmer_size, int32_t gdp_max_ed);
void printResult( gfa_edrst_t* r);

/*
 * This function writes results to a buffer. (Probably score  and maybe wlen.
 * Definitely path, which is the main output)
 * @param gfa_edrst_t r: the result of a kernel operation
 */
void storeResult(gfa_edrst_t* r);

/*
 * This function is called after all calls have been made to the kernel. It
 * invokes an expensive write operation writing all the files that you've stored
 * to buffers so far.
 * All other functions for storing store to local buffers in memory. When you
 * call this function it flushes the buffers to the filesystem.
 * This function must only be called once for the duration of the program. If it
 * is called more than once an error will be thrown
 */
void flushKernelDataBuffs();
/*
 * This function is used to store results produced by the kernel to the output
 * file, "
 *
 */
void storeResult(gfa_edrst_t* r);

/*
 * Sorry for the copy paste. Functions for storing scalar inputs to one line
 * each in a string buffer. 
 * V0 is the start vertex
 */
void storeV0(uint32_t v0);
/* V1 is end vertex */
void storeV1(uint32_t v1);
/* End0 is the end of the first linear chain (this is refference offset in node
 * I believe */
void storeEnd0(int32_t end0);
/* End1 is the start of the end linear chain (node offset I think) */
void storeEnd1(int32_t end1);
/* QueryGapLength is the distance in query between anchors you're bridging */
void storeQueryGapLength(uint32_t queryGapLen);
/* This stores the actual query sequence. The format is one query per line. Each
 * bp is one character, so we don't bother to add a delimiter between chars*/
void storeQueryGap(char* query, uint32_t gapLen);

#ifdef __cplusplus
}
#endif

#endif
