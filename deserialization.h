#ifndef DESERIALIZATION_H
#define DESERIALIZATION_H
/*
 * This file is intented to be very hacky. Not reusable. It's purpose is just a
 * quick sanity check to ensure I haven't damaged anything in the loading
 */
#include "gfa.h"

/*These are inverses of the serialization.h.*/
#ifdef __cplusplus
#include <cstdint>
extern "C" {
#endif
uint32_t loadV0();
uint32_t loadV1();
int32_t  loadEnd0();
int32_t  loadEnd1();
uint32_t loadQueryGapLen();
char*    loadQueryGap();
void initLoad();
gfa_t* loadGfa();
//This is direct copy from loadParams, but has C wrapper
const gfa_edseq_t* loadEsT(gfa_t* g);
void* loadKm();
#ifdef __cplusplus
}
#endif

#endif //DESERIALIZATION_H
