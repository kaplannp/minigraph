#ifndef LOAD_PARAMS_H
#define LOAD_PARAMS_H
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include "gfa.h"
#include "kalloc.h"

const int32_t GDP_MAX_ED = 10000;
//const KMER_SIZE = 19;

/*
 * Generic function for loading inputs if the format of the file is one input
 * per line of a simple type (e.g. uint32_t or int32_t.)
 */
template<typename T>
std::vector<T>* loadScalarInput(std::string filepath){
  std::vector<T>* inputs = new std::vector<T>();
  std::string line("");
  std::ifstream f(filepath);
  while (std::getline(f, line)){
    inputs->push_back(static_cast<T>(std::stoll(line)));
  }
  return inputs;
}

/*
 * The previous function doesn't quite generalize for strings, so we build
 * essentially a copy here that will do strings
 */
std::vector<std::string>* loadStrInput(std::string filepath);

/*
 * initializes es from the graph
 */
const gfa_edseq_t* loadEs(gfa_t* g);

/*
 * initialize a bunch of kalloc memory managers. I don't mean to be opaque, but
 * I have little inclination as to what these km objects actually do.
 * @param size the number of km objects to initialize
 */
std::vector<void*>* initKms(uint32_t size);

#endif//  LOAD_PARAMS_H
