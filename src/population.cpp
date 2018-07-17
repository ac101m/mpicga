// Standard headers
#include "unistd.h"
#include <iostream>
using namespace std;


// Project headers
#include "mpicga.hpp"
#include "mpi.h"
#include "utils.hpp"



// Construct the population with default constructors
population::population(uint32_t subPopulationCount, uint32_t genomeCount, uint32_t genomeLength) {

    // Initialise the algorithm
    this->algorithm = populationAlgorithm(subPopulationCount, genomeCount, genomeLength);

    // Population does not start initialised
    this->initialised = false;
}



// Initialise a population
void population::initialise(genomeTarget& target, uint32_t(*ff)(genomePerf_t)) {

    // Do population initialisation
    for(unsigned i = 0; i < this->algorithm.getSubPopulationCount(); i ++) {

        // Create the subpopulation and seed the internal random number generator
        this->subPopulations.push_back(subPopulation(this->algorithm.getSubPopulationAlgorithm()));
        this->subPopulations[i].getAlgorithm().setSeed(this->algorithm.localRand(0, (1 << 30) - 1));

        // Initialise the subpopulation
        this->subPopulations[i].initialise(target, ff, i);
    }

    // Initialise the sub population count vector
    // Vector contains the count of subPopulations resident on each rank
    for(int i = 0; i < rankCount(); i++) {

        // Calculate local subpopulation count
        uint32_t subPopulationCount = 0;
        for(unsigned j = 0; j < this->subPopulations.size(); j++) {
            if(domainDecomposition(j) == i) {
                subPopulationCount++;
            }
        }

        // Add the value of the counter to the subPopulation count vector
        this->rankSubPopulationCounts.push_back(subPopulationCount);
    }

    // Build an initial rankmap set initial fitness equal to
    // the rankmap index. This allows stability of sorting
    for(unsigned i = 0; i < this->subPopulations.size(); i++) {
        this->rankMap.push_back({&this->subPopulations[i], i, 0});
    }

    // Sort the initial rankmap
    this->updateRankMap();

    // We are now initialised
    this->initialised = true;
}



// Errors out if the population is not initialised
void population::assertInitialised(string msg) {
    if(!this->initialised) {
        errorOut(msg);
    }
}



// Gets pointers to all local subpopulations
vector<uint32_t> population::getLocalSubPopulationIndices(void) {

    // List of local virtual sub populations
    vector<uint32_t> localSubPopulationIndices;

    // Iterate over the array of
    for(unsigned i = 0; i < this->subPopulations.size(); i++) {
        if(this->subPopulations[i].isLocal()) {
            localSubPopulationIndices.push_back(i);
        }
    }

    // Return the list of local subPopulation
    return localSubPopulationIndices;
}



// Iterate the population n generations
void population::iterateSubPopulations(genomeTarget& target, uint32_t(*ff)(genomePerf_t), uint32_t n) {

    // Get a list of all local sub populations
    vector<uint32_t> localSubPopulationIndices = this->getLocalSubPopulationIndices();

    // Iterate all local subpopulations
    #pragma omp parallel for num_threads(12)
    for(unsigned i = 0; i < localSubPopulationIndices.size(); i++) {
        this->subPopulations[localSubPopulationIndices[i]].iterate(target, ff, n);
    }
}



// Swaps two indices in the fitness map
void population::swapRankMap(uint32_t i1, uint32_t i2) {

    // Temporary store
    auto tmp = this->rankMap[i1];

    // Perform the swap
    this->rankMap[i1] = this->rankMap[i2];
    this->rankMap[i2] = tmp;
}



// Get key to sort rank map member by
int64_t population::rankMapSortKey(uint32_t i) {
    return ((uint64_t)rankMap[i].fitness << 32) + rankMap[i].ptr->getDomainIndex();
}



// Partition function for the quicksort
int32_t population::partitionRankMap(int32_t low, int32_t high) {
    int64_t pivot = this->rankMapSortKey(high);
    int32_t i = low - 1;
    for(int32_t j = low; j <= high - 1; j++) {
        if(this->rankMapSortKey(j) <= pivot) {
            i++;
            this->swapRankMap(j, i);
        }
    }
    this->swapRankMap(i + 1, high);
    return i + 1;
}



// Perform the quicksort
void population::quickSortRankMap(int32_t low, int32_t high) {
    if(low < high) {
        int32_t pi = this->partitionRankMap(low, high);
        this->quickSortRankMap(low, pi - 1);
        this->quickSortRankMap(pi + 1, high);
    }
}



// Sorts the local copy of the rankmap
void population::sortRankMap(void) {
    quickSortRankMap(0, this->rankMap.size() - 1);
}



// Build a vector of integers respresentative of the local subpopulation
// Format: index - fitness - index - fitness - etc etc
// Remember to deallocate this when you are done!
unsigned *population::getRankMapTxBuffer(void) {

    // Create tx buffer structure
    unsigned *txBuffer = new unsigned[2 * this->getLocalSubPopulationCount()];

    // Get local subpopulation indices
    vector<unsigned> localSubPopIndices = this->getLocalSubPopulationIndices();

    // Iterate over local subpopulations and their indices to the buffer
    for(unsigned i = 0; i < this->getLocalSubPopulationCount(); i++) {

        // Calculate fitness (placeholder here)
        uint32_t fitness = this->subPopulations[localSubPopIndices[i]].getPerfData().bestGenomeFitness;

        // Push back index and fitness
        txBuffer[i * 2] = localSubPopIndices[i];
        txBuffer[(i * 2) + 1] = fitness;
    }

    // Return
    return txBuffer;
}



// Parses the rank map recieve buffer
void population::parseRankMapRxBuffer(unsigned *rxBuffer) {

    // Iterate over rankmap and fill from the recieve buffer
    for(unsigned i = 0; i < this->rankMap.size(); i++) {

        // Formulate rank map entry structure
        subPopulationFitnessMapping_t rankMapEntry;
        rankMapEntry.ptr = &this->subPopulations[rxBuffer[i * 2]];
        rankMapEntry.fitness = rxBuffer[(i * 2) + 1];

        // Add the rankmap entry to the
        this->rankMap[i] = rankMapEntry;
    }
}



// Synchronises the rankmap across all running processes
// I really need to clean this up
void population::synchroniseRankMap(void) {

    // Get local subpopulation indices
    unsigned *txBuffer = this->getRankMapTxBuffer();

    // Allocate raw memory for the recieve buffer
    unsigned *rxBuffer = new unsigned[2 * this->algorithm.getSubPopulationCount()];

    // Recieve buffer offsets
    vector<int> rxCounts;
    for(int i = 0; i < rankCount(); i++) {
        rxCounts.push_back(this->getSubPopulationCount(i) * 2);
    }

    // Build list of rank recieve counts
    vector<int> rxOffsets;
    unsigned currentOffset = 0;
    for(int i = 0; i < rankCount(); i++) {
        rxOffsets.push_back(currentOffset);
        currentOffset += this->rankSubPopulationCounts[i] * 2;
    }

    // Fill the recieve buffer with MPI_allgather
    MPI_Allgatherv(&txBuffer[0],
                   this->getLocalSubPopulationCount() * 2,
                   MPI_UNSIGNED,
                   &rxBuffer[0],
                   &rxCounts[0],
                   &rxOffsets[0],
                   MPI_UNSIGNED,
                   MPI_COMM_WORLD);

    // Parse the recieve buffer
    this->parseRankMapRxBuffer(rxBuffer);

    // Free the transimt and recieve buffers
    delete [] txBuffer;
    delete [] rxBuffer;
}



// Updates the rankmap
void population::updateRankMap(void) {

    // If multiple ranks are present, synchronise the rankmap across the processes
    this->synchroniseRankMap();

    // Sort the local copy of the rankmap
    this->sortRankMap();
}



// Get local subpopulation count
uint32_t population::getLocalSubPopulationCount(void) {
    return this->getSubPopulationCount(myRank());
}



// Get sub population count
uint32_t population::getSubPopulationCount(uint32_t rankAddress) {
    return this->rankSubPopulationCounts[rankAddress];
}



// Perform single crossover event
void population::doSubPopulationCrossover(genomeTarget& target, uint32_t(*ff)(genomePerf_t)) {

    // Do this the apropriate number of times
    for(unsigned i = 0; i < this->algorithm.getSelectCount(); i++) {

        // Generate population indices
        uint32_t pop1Idx = this->algorithm.randomHighSubPopulation();
        uint32_t pop2Idx = this->algorithm.randomHighSubPopulation();
        uint32_t destIdx = this->algorithm.randomLowSubPopulation();

        // Select two high and one low population
        subPopulation& pop1 = *this->rankMap[pop1Idx].ptr;
        subPopulation& pop2 = *this->rankMap[pop2Idx].ptr;
        subPopulation& destPop = *this->rankMap[destIdx].ptr;

        // Perform crossover
        destPop.crossover(pop1, pop2, this->algorithm.randomCrossoverIndices());
        destPop.updateRankMap(target, ff);

        // MPI barrier, hopefuly temporary
        MPI_Barrier(MPI_COMM_WORLD);
    }
}



// Iterate the population through one cycle
void population::iterate(genomeTarget& target, uint32_t(*ff)(genomePerf_t)) {

    // Make sure the population is initialised
    this->assertInitialised("Error, attempted to iterate uninitialised population.");

    // Do subpopulation crossover
    this->doSubPopulationCrossover(target, ff);

    // Iterate all local subpopulations by the apropriate number of generations per cycle
    this->iterateSubPopulations(target, ff, this->algorithm.getGenerationsPerCycle());

    // Synchonise the global rankmap across all processes
    this->updateRankMap();
}



// Iterate the population through n cycles
void population::iterate(genomeTarget& target, uint32_t(*ff)(genomePerf_t), uint32_t n) {

    // Iterate the population through n cycles
    for(unsigned i = 0; i < n; i++) {
        this->iterate(target, ff);
        this->rankMap[0].ptr->printRankMap(target);
    }
}



// Print the rank map on each of the ranks indicated in "ranks"
void population::printRankMap(void) {

    // Iterate over all processes
    for(int i = 0; i < rankCount(); i++) {

        // Serialise ranks
        if(i == myRank()) {

            // Print out rank name
            cout << rankString() << "\n";

            // Print out the rank map
            for(unsigned j = 0; j < this->rankMap.size(); j++) {
                cout << "Ranking: " << j;
                cout << " Index: " << this->rankMap[j].ptr->getDomainIndex();
                cout << " Fitness: " << this->rankMap[j].fitness;
                cout << "\n";
            }
            cout << "\n";
        }

        // Wait for previous rank to finish first
        usleep(20000);
        MPI_Barrier(MPI_COMM_WORLD);
    }
}
