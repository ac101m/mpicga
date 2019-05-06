#ifndef UTILS_H
#define UTILS_H


// Standard
#include "stdint.h"
#include <string>


// Counts bits in a 64 bit word
uint32_t countBits(uint64_t data);


// Generates rank std::string for this rank
std::string rankString(void);


// Oh shit son, something done fucked up
void err(std::string msg);


// Ehh, something is a bit squiffy
void warn(std::string msg);


// Get the address of this rank
int32_t myRank(void);


// Get the number of rancs in the communicator
int32_t rankCount(void);



#endif // UTILS_H
