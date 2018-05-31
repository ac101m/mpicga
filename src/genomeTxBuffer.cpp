// Standard headers
#include <iostream>


// Project headers
#include "mpicga.hpp"
#include "utils.hpp"



// Constructor
genomeTransmissionBuffer::genomeTransmissionBuffer(unsigned bufferLength) {
    this->buffer = new geneNetworkFrame_t[bufferLength];
    this->maxGenes = bufferLength;
    this->currentGenes = 0;
}



// Get the raw data in the form of a vector of gene network frames
geneNetworkFrame_t *genomeTransmissionBuffer::getData(void) {

    // Return the gene data
    return &this->buffer[0];
}



// Append a gene network frame
void genomeTransmissionBuffer::append(geneNetworkFrame_t g) {

    // Check that buffer is not full
    if(currentGenes == maxGenes) {
        errorOut("Error, genome transmit buffer overflow (append).");
    }

    // Position in the buffer to start writing at
    this->buffer[this->currentGenes] = g;

    // Increment the gene count
    this->currentGenes++;
}



// Append a gene
void genomeTransmissionBuffer::append(gene g) {

    // Add it to the internal vector of gene network frames
    this->append(g.getNetworkFrame());
}



// Append a genome
void genomeTransmissionBuffer::append(genome g) {

    // Add it to the internal vector of gene network frames
    for(unsigned i = 0; i < g.getGeneCount(); i++) {
        this->append(g.getGenes()[i]);
    }
}



// Transmit the buffer
void genomeTransmissionBuffer::transmit(int32_t destination, int32_t tag) {

    // Transmit the buffer using MPI
    MPI_Ssend((uint8_t *)&this->buffer[0],
             currentGenes * sizeof(geneNetworkFrame_t),
             MPI_BYTE,
             destination,
             tag,
             MPI_COMM_WORLD);
}



// Recieve buffer from another transmit
void genomeTransmissionBuffer::receive(int32_t source, int32_t tag) {

    // Probe and get status of recieve operation
    MPI_Status status;
    MPI_Probe(source,
              tag,
              MPI_COMM_WORLD,
              &status);

    // Get total number of bytes that are to be transmitted
    int32_t byteCount;
    MPI_Get_count(&status, MPI_BYTE, &byteCount);

    // Confirm that incoming buffer is not malformed
    if(byteCount % sizeof(geneNetworkFrame_t) != 0) {
        errorOut("Error, recieve into buffer of incorrect size.");
    } else {
        if(this->maxGenes != (byteCount / sizeof(geneNetworkFrame_t))) {
            errorOut("Error, genome recieve buffer is of incorrect length.");
        } else {
            this->currentGenes = byteCount / sizeof(geneNetworkFrame_t);
        }
    }

    // Recieve the buffer using MPI
    MPI_Recv((uint8_t *)&this->buffer[0],
             byteCount,
             MPI_BYTE,
             source,
             tag,
             MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
}


// Destructor
genomeTransmissionBuffer::~genomeTransmissionBuffer(void) {
    delete [] this->buffer;
}
