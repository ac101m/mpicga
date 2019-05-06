// Standard headers
#include <iostream>
#include <vector>
using namespace std;


// Project headers
#include "utils.hpp"
#include "mpicga.hpp"



// Constructor, initialises everything to zero
gene::gene(void) {
  this->bufValid = false;
  this->function = GENE_FN_NOP;
  this->aIndex = 0;
  this->bIndex = 0;
  this->buf = 0;
}



// Construct a gene from a gene network frame
gene::gene(geneNetworkFrame_t frame) {
  this->bufValid = false;
  this->function = frame.function;
  this->aIndex = frame.aIndex;
  this->bIndex = frame.bIndex;
  this->buf = 0;
}



// Gets output buffer, recomputes if neccessary
uint64_t gene::getOutputBuffer(vector<gene>& genes) {
  uint64_t aInput, bInput = 0;

  // Check that the input buffer is valid
  if(!this->bufValid) {

    // Recursively evaluate gene values
    aInput = genes[this->aIndex].getOutputBuffer(genes);
    if((this->function != GENE_FN_NOP) && (this->function != GENE_FN_NOT)) {
      bInput = genes[this->bIndex].getOutputBuffer(genes);
    }

    // Mark gene output as valid
    this->buf = this->computeBufferValue(aInput, bInput);
    this->bufValid = true;
  }

  // Return the output buffer
  return this->buf;
}



// Compute buffer value
uint64_t gene::computeBufferValue(uint64_t a, uint64_t b) {

  // Compute value depending on
  switch(this->function) {

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
      err("Error, encountered unrecognised gene function during evaluation.\n");
      return 0;
      break;
  }
}



// Function to set gene input
void gene::overrideBuffer(uint64_t bv) {
  this->buf = bv;
  this->bufValid = true;
}



// Function to randomly mutate the gene
bool gene::mutate(uint32_t selectedIndex, subPopulationAlgorithm& algorithm) {

  // Randomly select gene characteristic to mutate
  switch(algorithm.localRand(0, 2)) {
    case 0: this->aIndex = algorithm.randomGeneInputIndex(selectedIndex); break;
    case 1: this->bIndex = algorithm.randomGeneInputIndex(selectedIndex); break;
    case 2: this->function = algorithm.randomGeneFunction(); break;
    default:
      err("Error, failed gene mutation operation.\n");
      exit(1);
      break;
  }

  // Invalidate output buffer
  bool previouslyActive = this->bufValid;
  this->bufValid = false;
  return previouslyActive;
}



// Get transmission template
geneNetworkFrame_t gene::getNetworkFrame(void) {
  geneNetworkFrame_t t;

  // Populate the gene template struct
  t.aIndex = this->aIndex;
  t.bIndex = this->bIndex;
  t.function = this->function;

  // Return the gene template struct
  return t;
}
