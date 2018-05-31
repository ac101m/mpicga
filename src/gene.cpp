// Standard headers
#include <iostream>
#include <vector>
using namespace std;


// Project headers
#include "utils.hpp"
#include "mpicga.hpp"



// Constructor, initialises everything to zero
gene::gene(void) {
    this->outputBufferValid = false;
    this->geneFunction = GENE_FN_NOP;
    this->aInputIndex = 0;
    this->bInputIndex = 0;
    this->outputBuffer = 0;
}



// Construct a gene from a gene network frame
gene::gene(geneNetworkFrame_t frame) {
    this->outputBufferValid = false;
    this->geneFunction = frame.geneFunction;
    this->aInputIndex = frame.aInputIndex;
    this->bInputIndex = frame.bInputIndex;
    this->outputBuffer = 0;
}



// Gets output buffer, recomputes if neccessary
uint64_t gene::getOutputBuffer(vector<gene>& genes) {
    uint64_t aInput, bInput = 0;

    // Check that the input buffer is valid
    if(!this->outputBufferValid) {

        // Recursively evaluate gene values
        aInput = genes[this->aInputIndex].getOutputBuffer(genes);
        if((this->geneFunction != GENE_FN_NOP) && (this->geneFunction != GENE_FN_NOT)) {
            bInput = genes[this->bInputIndex].getOutputBuffer(genes);
        }

        // Mark gene output as valid
        this->outputBuffer = this->computeBufferValue(aInput, bInput);
        this->outputBufferValid = true;
    }

    // Return the output buffer
    return this->outputBuffer;
}



// Compute buffer value
uint64_t gene::computeBufferValue(uint64_t a, uint64_t b) {

    // Compute value depending on
    switch(this->geneFunction) {

        // Logic functions
        case GENE_FN_NOP: return a; break;
        case GENE_FN_NOT: return ~a; break;
        case GENE_FN_AND: return a & b; break;
        case GENE_FN_NAND: return ~(a & b); break;
        case GENE_FN_OR: return a | b; break;
        case GENE_FN_NOR: return ~(a | b); break;
        case GENE_FN_XOR: return a ^ b; break;
        case GENE_FN_XNOR: return ~(a ^ b); break;

        // Unrecognised function
        default:
            errorOut("Error, encountered unrecognised gene function during evaluation.\n");
            return 0;
            break;
    }
}



// Function to set gene input
void gene::overrideBuffer(uint64_t bv) {
    this->outputBuffer = bv;
    this->outputBufferValid = true;
}



// Function to randomly mutate the gene
bool gene::mutate(uint32_t selectedIndex, subPopulationAlgorithm& algorithm) {

    // Randomly select gene characteristic to mutate
    switch(algorithm.localRand(0, 2)) {
        case 0: this->aInputIndex = algorithm.randomGeneInputIndex(selectedIndex); break;
        case 1: this->bInputIndex = algorithm.randomGeneInputIndex(selectedIndex); break;
        case 2: this->geneFunction = algorithm.randomGeneFunction(); break;
        default:
            errorOut("Error, failed gene mutation operation.\n");
            exit(1);
            break;
    }

    // Invalidate output buffer
    bool previouslyActive = this->outputBufferValid;
    this->outputBufferValid = false;
    return previouslyActive;
}



// Get transmission template
geneNetworkFrame_t gene::getNetworkFrame(void) {
    geneNetworkFrame_t t;

    // Populate the gene template struct
    t.aInputIndex = this->aInputIndex;
    t.bInputIndex = this->bInputIndex;
    t.geneFunction = this->geneFunction;

    // Return the gene template struct
    return t;
}
