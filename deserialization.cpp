#include "loadParams.h"
#include "deserialization.h"
#include <vector>
#include "gfa.h"

std::vector<uint32_t>* V0Buff;
std::vector<std::string>* QueryGapBuff;
std::vector<uint32_t>* QueryGapLenBuff;
std::vector<int32_t>* End1Buff;
std::vector<int32_t>* End0Buff;
std::vector<uint32_t>* V1Buff;
std::vector<void*>* KmBuff;
gfa_t *gfa;

void initLoad(){
  V0Buff = loadScalarInput<uint32_t>("DumpToLoad/Inputs/v0.txt");
  QueryGapBuff = loadStrInput("DumpToLoad/Inputs/queryGap.txt");
  QueryGapLenBuff = loadScalarInput<uint32_t>("DumpToLoad/Inputs/queryGapLen.txt");
  End1Buff = loadScalarInput<int32_t>("DumpToLoad/Inputs/end1.txt");
  End0Buff = loadScalarInput<int32_t>("DumpToLoad/Inputs/end0.txt");
  V1Buff = loadScalarInput<uint32_t>("DumpToLoad/Inputs/v1.txt");
  gfa = gfa_read("test/MT.gfa");
  KmBuff = initKms(V0Buff->size());
}

uint32_t loadV0(){
  static int i=0;
  return (*V0Buff)[i++];
}
uint32_t loadV1(){
  static int i=0;
  return (*V1Buff)[i++];
}
int32_t loadEnd0(){
  static int i=0;
  return (*End0Buff)[i++];
}
int32_t  loadEnd1(){
  static int i=0;
  return (*End1Buff)[i++];
}
uint32_t loadQueryGapLen(){
  static int i=0;
  return (*QueryGapLenBuff)[i++];
}
char*    loadQueryGap(){
  static int i=0;
  return (*QueryGapBuff)[i++].data();
}
gfa_t* loadGfa(){
  return gfa;
}

const gfa_edseq_t* loadEsT(gfa_t* g){
  loadEs(g); //call from loadParams
}

void* loadKm(){
  static int i=0;
  return (*KmBuff)[i++];
}
