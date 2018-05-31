#ifndef UTILS_H
#define UTILS_H


// Standard c stuff
#include "stdint.h"
#include <string>
using namespace std;


// Counts bits in a 64 bit word
uint32_t countBits(uint64_t data);


// Generates rank string for this rank
string rankString(void);


// Oh shit son, something done fucked up
void errorOut(string msg);


// Ehh, something is a bit squiffy
void warningOut(string msg);


// Get the address of this rank
int32_t myRank(void);


// Get the number of rancs in the communicator
int32_t rankCount(void);



#endif // UTILS_H
