#include "loadParams.h"
#include "gfa.h"

std::vector<std::string>* loadStrInput(std::string filepath){
  std::vector<std::string>* inputs = new std::vector<std::string>();
  std::string line("");
  std::ifstream f(filepath);
  while (std::getline(f, line)){
    inputs->push_back(line);
  }
  return inputs;
}

const gfa_edseq_t* loadEs(gfa_t* g){
  return gfa_edseq_init(g);
}

std::vector<void*>* initKms(uint32_t size){
  std::vector<void*>* kms = new std::vector<void*>();
  for (uint32_t i = 0; i < size; i++){
    kms->push_back(km_init());
  }
  return kms;
}
