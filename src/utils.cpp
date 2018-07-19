// This sources header
#include "utils.hpp"


// CPP
#include <sstream>
#include <iostream>


// C
#include "stdlib.h"
#include "stdint.h"


// MPI
#include "mpi.h"



// Naive set bit counter, will optimise later (maybe)
uint32_t countBits(uint64_t data) {
    uint32_t bitCount = 0;

    uint64_t bits = data;
    while(bits) {
        bitCount++;
        bits &= (bits - 1);
    }

    return bitCount;
}



// Returns rank address of running process
int32_t myRank(void) {
    int32_t myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    return myRank;
}



// Returns number of
int32_t rankCount(void) {
    int32_t rankCount;
    MPI_Comm_size(MPI_COMM_WORLD, &rankCount);
    return rankCount;
}



// Generates rank string for this rank
string rankString(void) {

    // Get processor name
    int32_t length;
    char procName[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(&procName[0], &length);

    // Get rank number
    int32_t processRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);

    // Build stringstream
    stringstream ss;
    ss << "[" << procName << ":" << processRank << "]";
    return ss.str();
}



// Oh shit son, something done fucked up
void err(string msg) {

    // Get processor name & rank string
    cout << rankString() << " " << msg;

    // Make sure message ends with newline
    if(msg[msg.size() - 1] != '\n') {
        cout << "\n";
    }

    // El fin
    MPI_Abort(MPI_COMM_WORLD, 0);
    exit(1);
}



// Oh shit son, something done fucked up
void warn(string msg) {

    // Get processor name & rank string
    cout << rankString() << " " << msg;

    // Make sure message ends with newline
    if(msg[msg.size() - 1] != '\n') {
        cout << "\n";
    }
}
