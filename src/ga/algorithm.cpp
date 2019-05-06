// Standard headers
#include <iostream>
using namespace std;


// Project headers
#include "mpicga.hpp"
#include "utils.hpp"



//========[SUBPOPULATION ALGORITHM]==============================================================//

// Subpopulation algorithm class default constructor
subPopulationAlgorithm::subPopulationAlgorithm(uint32_t genomeCount, uint32_t genomeLength) {

  // Seed the random number generator
  this->setSeed(0);

  // Population geometry
  this->genomeCount = genomeCount;
  this->genomeLength = genomeLength;

  // Default allowable gene functions
  this->allowableFunctions = {GENE_FN_AND,
                              GENE_FN_OR,
                              GENE_FN_XOR,
                              GENE_FN_NOT};

  // Default selection and mutation counts
  this->mutateCount = 1;
  this->selectCount = 1;

  // High and low select ranges
  this->highSelectRange = this->lowSelectRange = this->genomeCount / 2;

  // Subpopulation min and max feed forward fractions
  this->minFeedForward = 1;
  this->maxFeedForward = this->getGenomeLength();
}



// Default constructor
subPopulationAlgorithm::subPopulationAlgorithm(void) {

  // Call the sized constructor
  subPopulationAlgorithm(8, 128);
}


// Get max gate delays
void subPopulationAlgorithm::setMinGateDelays(uint32_t gd) {

  // If min gate delays is set to zero, set max feed forward
  // to match the length of the genome
  if(!gd) {
    this->maxFeedForward = this->genomeLength;
  }

  // Else, set up accordingly
  this->maxFeedForward = this->genomeLength / gd;
}



// Get min gate delays
void subPopulationAlgorithm::setMaxGateDelays(uint32_t gd) {
  this->minFeedForward = this->genomeLength / gd;
}



// Rnadom number generation using local number generator
int32_t subPopulationAlgorithm::localRand(int32_t minimum, int32_t maximum) {

  // Better random code
  uniform_int_distribution<> distribution(minimum, maximum);
  return distribution(this->localRandEngine);
}



// Generate random high genome index
// Fit genomes at low indexes (0 = most fit)
int32_t subPopulationAlgorithm::randomHighGenome(void) {
  int32_t rand = this->highSelectRange - 1;

  // Generate the random index
  for(unsigned i = 0; i < 2; i++) {
    rand = this->localRand(0, rand);
  }

  // Return the randomly generated index
  return rand;
}



// Generate random low genome index
// Unfit genomes at high indexes
int32_t subPopulationAlgorithm::randomLowGenome(void) {
  int32_t rand = this->lowSelectRange - 1;

  // Generate the random index
  for(unsigned i = 0; i < 2; i++) {
    rand = this->localRand(0, rand);
  }

  // Return the randomly generated index
  return (this->genomeCount - 1) - rand;
}



// Generate a random index apropriate to a gene at position i
uint32_t subPopulationAlgorithm::randomGeneInputIndex(int32_t i) {

  // Variables for start and end range ov valid indices
  int32_t rangeStart = i - this->getMaxFeedForward();
  int32_t rangeEnd = i - this->getMinFeedForward();

  // Make sure valid range does not start somewhere
  // "off the bottom" of the gene vector
  if(rangeStart < 0) {
    rangeEnd += (0 - rangeStart);
    rangeStart = 0;
  }

  // Make sure valid range does not exceed selected gene index
  // (graph represented by genome *must* remain acyclic)
  if(rangeEnd >= i) {
    rangeEnd = i - 1;
  }

  // Generate and return the index
  return this->localRand(rangeStart, rangeEnd);
}



// Pick a gene at random
geneFunction_t subPopulationAlgorithm::randomGeneFunction(void) {
  return this->allowableFunctions[this->localRand(0, this->allowableFunctions.size() - 1)];
}



//========[POPULATION ALGORITHM]=================================================================//

// Default constructor
populationAlgorithm::populationAlgorithm(void) {

  // Call the initialising constructor with default values
  populationAlgorithm(4, 8, 128);
}



// Initialising constructor
populationAlgorithm::populationAlgorithm(uint32_t subPopCount, uint32_t genomeCount, uint32_t genomeLength) {

  // Initialise the subpopulation algorithm and count
  this->subPopAlgorithm = subPopulationAlgorithm(genomeCount, genomeLength);

  // Initialise subpopulation count
  this->subPopulationCount = subPopCount;

  // Default generations per cycle
  this->generationsPerCycle = 65536;

  // Set the seed to 1 by default
  this->localRandEngine.seed(1);

  // Crossover related variables
  this->selectCount = 1;
  this->crossoverCount = 4;
  this->highSelectRange = this->lowSelectRange = this->subPopulationCount / 2;

  // Processing behaviour
  this->threadCount = 1;
  this->commTagCounter = 0;
}



// Random number generation using local number generator
int32_t populationAlgorithm::localRand(int32_t minimum, int32_t maximum) {

  // Better random code
  uniform_int_distribution<> distribution(minimum, maximum);
  return distribution(this->localRandEngine);
}



// Selects a low subPopulation index at random
int32_t populationAlgorithm::randomHighSubPopulation(void) {
  int32_t rand = this->highSelectRange - 1;

  // Iterate, for exponential selection
  for(unsigned i = 0; i < 2; i++) {
    rand = this->localRand(0, rand);
  }

  // Return the randomly generated number
  return rand;
}



// Selects a high subPopulation index at random
int32_t populationAlgorithm::randomLowSubPopulation(void) {
  int32_t rand = this->lowSelectRange - 1;

  // Iterate, for exponential selection
  for(unsigned i = 0; i < 2; i++) {
    rand = this->localRand(0, rand);
  }

  // Return the randomly generated number
  return (this->subPopulationCount - 1) - rand;
}



// Generates a random list of crossover indices
vector<uint32_t> populationAlgorithm::randomCrossoverIndices(void) {

  // Build vector of crossover points
  vector<uint32_t> indices;
  for(unsigned i = 0; i < this->crossoverCount; i++) {
    indices.push_back(this->localRand(0, this->subPopulationCount - 1));
  }

  // Return the crossover indices vector
  return indices;
}
