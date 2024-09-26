#include "serialization.h"
#include "minigraph.h"
#include <iostream>
#include <vector>
#include <cstdint>
#include <sstream>
#include <fstream>

//We keep a bunch of global variables here and they are all dumped to their
//respective files at the end
std::stringstream resultBuff;
std::stringstream v0Buff;
std::stringstream v1Buff;
std::stringstream end0Buff;
std::stringstream end1Buff;
std::stringstream queryGapLenBuff;
std::stringstream queryGapBuff;

void storeV0(uint32_t v0){
  v0Buff << v0 << std::endl;
}
void storeV1(uint32_t v1){
  v1Buff << v1 << std::endl;
}
void storeEnd0(int32_t end0){
  end0Buff << end0 << std::endl;
}
void storeEnd1(int32_t end1){
  end1Buff << end1 << std::endl;
}
void storeQueryGapLength(uint32_t queryGapLen){
  queryGapLenBuff << queryGapLen << std::endl;
}

void storeQueryGap(char* query, uint32_t queryGapLen){
  for(uint32_t i = 0; i < queryGapLen; i++){
    queryGapBuff << query[i];
  }
  queryGapBuff << std::endl;
}

void storeResult(gfa_edrst_t* r){
  //TODO. I would remove these first two as they just confuse, but want to be
  //absolutely sure that my version is correct first
  //push the s variable (which we think is score. if the author was consisten in
  //labelling)
  resultBuff << r->s << " ";
  //push the wlen variable which appears to be some kind of character count. I
  resultBuff << r->wlen << " ";

  resultBuff << "{ ";
  //std::cerr << r->nv << std::endl;
  for(int i = 0; i < r->nv; i++){
    resultBuff << r->v[i] << " ";
  }
  resultBuff << "}";
  resultBuff << std::endl;
}

/*
 * helper funtion for writeAll which deletes the directory `Dump` and creates
 *   Dump
 *   Dump/Inputs
 *   Dump/Out
 * and will exit program if one of these operation fails
 */
void initializeDirs(){
  const std::string dumpDir = "Dump";
  const std::string inputDir = dumpDir + "/Inputs";
  const std::string outDir = dumpDir + "/Out";
  bool failed = false;
  std::string stager("");
  //initialize the directories if this is the first time in the function
  stager = "rm -rf " + dumpDir;
  failed |= std::system(stager.c_str());
  stager = "mkdir " + dumpDir;
  failed |= std::system(stager.c_str());
  stager = "mkdir " + inputDir;
  failed |= std::system(stager.c_str());
  stager = "mkdir " + outDir;
  failed |= std::system(stager.c_str());
  if (failed){
    std::cerr << "initializing dump directory failed! Aborting" << std::endl;
    exit(1);
  }
}

void flushKernelDataBuffs(){
  static bool called = false;
  if (called){
    std::cerr << "writeAll called twice! This function must be called only once"
              << " at the end of the kernel";
    exit(1);
  }
  initializeDirs();
  const std::string dumpDir = "Dump";
  const std::string inputDir = dumpDir + "/Inputs";
  const std::string outDir = dumpDir + "/Out";
  std::ofstream resultsDir(outDir+"/results.txt");
  std::ofstream v0Dir(inputDir+"/v0.txt");
  std::ofstream v1Dir(inputDir+"/v1.txt");
  std::ofstream end0Dir(inputDir+"/end0.txt");
  std::ofstream end1Dir(inputDir+"/end1.txt");
  std::ofstream queryGapLenDir(inputDir+"/queryGapLen.txt");
  std::ofstream queryGapDir(inputDir+"/queryGap.txt");

  resultsDir << resultBuff.str();
  v0Dir << v0Buff.str();
  v1Dir << v1Buff.str();
  end0Dir << end0Buff.str();
  end1Dir << end1Buff.str();
  queryGapLenDir << queryGapLenBuff.str();
  queryGapDir << queryGapBuff.str();
  called=true;
}

void print(
    bridge_aux_t* aux,
    int32_t kmer_size,
    int32_t gdp_max_ed
    ){

  std::cerr << " aux.n_llc=" << aux->n_llc;
  std::cerr << " aux.m_llc=" << aux->m_llc;
  std::cerr << " aux.n_a=" << aux->n_a;
  std::cerr << " aux.n_seg=" << aux->n_seg;
  std::cerr << " aux.es=";
  for(int i=0; i < aux->es->len; i++){
    std::cerr << aux->es->seq[i];
  }
  std::cerr << " kmer_size=" << kmer_size;
  std::cerr << " gdp_max_ed=" << gdp_max_ed << std::endl;

}

void printResult( gfa_edrst_t* r){
  std::cerr << " s=" << r->s;
  std::cerr << " end_v=" << r->end_v;
  std::cerr << " end_off=" << r->end_off;
  std::cerr << " wlen=" << r->wlen;
  std::cerr << " n_end=" << r->n_end;
  std::cerr << " nv=" << r->nv;
  std::cerr << " n_iter=" << r->n_iter;

  std::cerr << " path={";
  int length=0;
  uint32_t node;
  do {
    node = r->v[length];
    std::cerr << node << ", ";
    length++; 
  } while (node != r->end_v);
  std::cerr << "} length=" << length << std::endl;
}



//void stKmer(int32_t 
