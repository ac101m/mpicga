// Standard headers
#include <iostream>
using namespace std;


// Project headers
#include "mpicga.hpp"
#include "utils.hpp"



// Population domain decomposition function
int32_t domainDecomposition(uint32_t indexWithinDomain) {

    // Get the total number of ranks
    int32_t rankCount;
    MPI_Comm_size(MPI_COMM_WORLD, &rankCount);

    // Calculate location rank
    return indexWithinDomain % rankCount;
}



// Constructor, builds with specific population size, but default algorithm
subPopulation::subPopulation(uint32_t populationSize, uint32_t genomeSize) {

    // Reset the subpopulation algorithm
    this->algorithm = subPopulationAlgorithm(populationSize, genomeSize);

    // Clear all data structures
    this->genomes.clear();
    this->rankMap.clear();

    // Initialise world address to zero
    this->commWorldAddress = 0;
    this->local = false;

    // This subpopulation is not initialised
    this->initialised = false;
}



// Constructor, builds with passed-in algorithm
subPopulation::subPopulation(subPopulationAlgorithm algorithm) {

    // Initialise the input algorithm
    this->algorithm = algorithm;

    // Clear all data structures
    this->genomes.clear();
    this->rankMap.clear();

    // Initialise world address to zero
    this->commWorldAddress = 0;
    this->local = false;

    // This subpopulation is not initialised
    this->initialised = false;
}



// Initialisation method
// Initialises subpopulation on indicated process
void subPopulation::initialise(genomeTarget& target, uint32_t(*ff)(genomePerf_t), int32_t domainIndex) {

    // Initialise comm world address
    this->domainIndex = domainIndex;
    this->commWorldAddress = domainDecomposition(this->domainIndex);

    // Is this genome local or not?
    if(this->commWorldAddress == myRank()) {

        // Initialise random genomes to the genome vector
        for(unsigned i = 0; i < algorithm.getGenomeCount(); i++) {
            this->genomes.push_back(genome(algorithm.getGenomeLength(), algorithm));
        }

        // Build the initial rankmap
        for(unsigned i = 0; i < this->genomes.size(); i++) {
            this->rankMap.push_back({&genomes[i], i, 0});
        }

        // Update the rankmap to contain fitness values for randomly generated genomes
        this->updateRankMap(target, ff);

        // Indicate that this subpopulation is local to this process
        this->local = true;

    } else {
        this->local = false;
    }

    // Subpopulation is now initialised
    this->initialised = true;
}



// Initialisation method
// Initialises on process zero
void subPopulation::initialise(genomeTarget& target, uint32_t(*ff)(genomePerf_t)) {

    // Initialise as if on rank 0
    this->initialise(target, ff, 0);
}



// Returns true if subpopulation is local to this node
bool subPopulation::isLocal(void) {

    // Assert that the population is initialised
    this->assertInitialised("Error, attempt to query locality of uninitialised subpopulation.");

    // Return locality of genome
    return this->local;
}



// Swaps two indices in the fitness map
void subPopulation::swapRankMap(uint32_t i1, uint32_t i2) {

    // Temporary store for a tuple
    auto tmp = this->rankMap[i1];

    // Perform the swap
    this->rankMap[i1] = this->rankMap[i2];
    this->rankMap[i2] = tmp;
}



// Get key to sort rank map with
int64_t subPopulation::rankMapSortKey(uint32_t i) {
    return ((uint64_t)rankMap[i].fitness << 32) + rankMap[i].index;
}



// Partition function for the quicksort
int32_t subPopulation::partitionRankMap(int32_t low, int32_t high) {
    int64_t pivot = this->rankMapSortKey(high);
    int32_t i = low - 1;
    for(int32_t j = low; j <= high - 1; j++) {
        if((int64_t)this->rankMapSortKey(j) <= pivot) {
            i++;
            this->swapRankMap(i, j);
        }
    }
    this->swapRankMap(i + 1, high);
    return i + 1;
}



// Quicksort, somehow
void subPopulation::quickSortRankMap(int32_t low, int32_t high) {
    if(low < high) {
        int32_t pi = this->partitionRankMap(low, high);
        this->quickSortRankMap(low, pi - 1);
        this->quickSortRankMap(pi + 1, high);
    }
}



// Sorts the things
void subPopulation::sortRankMap(void) {
    this->quickSortRankMap(0, this->rankMap.size() - 1);
}



// Updates the rankmap
void subPopulation::updateRankMap(genomeTarget& target, uint32_t(*ff)(genomePerf_t)) {

    // Update the rankmap fitness values
    for(unsigned i = 0; i < this->rankMap.size(); i++) {
        this->rankMap[i].fitness = ff(this->rankMap[i].ptr->getPerfData(target));
    }

    // Sort the rankmap
    this->sortRankMap();
}



// Iterate the population using specific mutation specs
void subPopulation::iterate(genomeTarget& target, uint32_t(*ff)(genomePerf_t)) {

    // Assert that the population is initialised
    this->assertInitialised("Error, attempted to iterate uninitialised subpopulation.");

    // For every selection
    for(unsigned i = 0; i < this->algorithm.getSelectCount(); i++) {

        // Select genome indices within the rankmap
        uint32_t fitIdx = this->algorithm.randomHighGenome();
        uint32_t unfitIdx = this->algorithm.randomLowGenome();

        // Get pointers to genomes
        genome* fitGenome = this->rankMap[fitIdx].ptr;
        genome* unfitGenome = this->rankMap[unfitIdx].ptr;

        // Select a genome, and mutate
        if (fitIdx != unfitIdx) {
            *unfitGenome = *fitGenome;
            unfitGenome->mutate(this->algorithm);
        }
    }

    // Increment genome ages
    for(unsigned i = 0; i < this->rankMap.size(); i++) {
        this->rankMap[i].ptr->incrementAge();
    }

    // Update the rankmap
    this->updateRankMap(target, ff);
}



// Iterates the population n times
void subPopulation::iterate(genomeTarget& target, uint32_t(*ff)(genomePerf_t), uint32_t n) {

    // Iterate the population n times
    for(unsigned i = 0; i < n; i++) {
        this->iterate(target, ff);
    }
}



// Throws an error if the subpopulation is not initialised
void subPopulation::assertInitialised(string msg) {
    if(!this->initialised) {
        errorOut(msg);
    }
}



// Throws an error of the calling rank is not the subpopulation local rank
void subPopulation::assertLocal(string msg) {
    if(!this->local) {
        errorOut(msg);
    }
}



// Get performance data for this subpopulation
subPopulationPerf_t subPopulation::getPerfData(void) {

    // Assert that things are going alright
    this->assertInitialised("Error, attempted to get performance data of uninitialised subpopulation.");
    this->assertLocal("Error, attempted to get performance data from nonlocal subpopulation.");

    // Populate the struct and return
    subPopulationPerf_t perf;
    perf.bestGenomeFitness = this->rankMap[0].fitness;

    // Return the performance data
    return perf;
}



// Copy all specified genome indices to
void subPopulation::copyGenomes(vector<uint32_t>& genomeIndices, subPopulation& source) {
    for(unsigned i = 0; i < genomeIndices.size(); i++) {
        this->genomes[genomeIndices[i]].copyFrom(source.getGenomes()[genomeIndices[i]]);
    }
}



// Parse buffer
void subPopulation::parseGenomeBuffer(genomeTransmissionBuffer& buffer, vector<uint32_t>& genomeIndices) {

    // Get the raw buffer
    geneNetworkFrame_t *geneData = buffer.getData();
    geneNetworkFrame_t *geneNetworkFramePtr;

    // Parse the genomes one by one
    for(unsigned i = 0; i < genomeIndices.size(); i++) {
        geneNetworkFramePtr = &geneData[i * this->algorithm.getGenomeLength()];
        this->genomes[genomeIndices[i]].parseGeneNetworkFrameArray(geneNetworkFramePtr);
    }
}



// Transmits genomes
void subPopulation::exportGenomes(vector<uint32_t>& genomeIndices, subPopulation& target) {

    // Make sure we are local
    this->assertLocal("Error, attempt to export genomes from nonlocal subpopulation.");

    // Create a transmit buffer
    genomeTransmissionBuffer txBuffer(genomeIndices.size() * this->algorithm.getGenomeLength());

    // Add genomes to export to the buffer
    for(unsigned i = 0; i < genomeIndices.size(); i++) {
        txBuffer.append(this->genomes[genomeIndices[i]]);
    }

    // Transmit the contents of the buffer
    txBuffer.transmit(target.getProcessRank(), this->getDomainIndex());
}



// Recieves genomes
void subPopulation::importGenomes(vector<uint32_t>& genomeIndices, subPopulation& source) {

    // Make sure we are local
    this->assertLocal("Error, attempt to import genomes to nonlocal subpopulation.");

    // Create a recieve buffer of apropriate size
    genomeTransmissionBuffer rxBuffer(source.getAlgorithm().getGenomeLength() * genomeIndices.size());

    // Perform the recieve operation
    rxBuffer.receive(source.getProcessRank(), source.getDomainIndex());

    // Parse the genomes from the input buffer
    this->parseGenomeBuffer(rxBuffer, genomeIndices);
}



// Subpopulation crossover operator
void subPopulation::crossover(subPopulation& pop1, subPopulation& pop2, vector<uint32_t> crossoverIndices) {

    // Check that this is initialised
    assertInitialised("Error, attempted to perform crossover operation on uninitialised supopulation.");

    // Vectors to contain vectors of transmit and recieve genome indices
    vector<uint32_t> p1Indices; p1Indices.reserve(this->algorithm.getGenomeCount());
    vector<uint32_t> p2Indices; p2Indices.reserve(this->algorithm.getGenomeCount());

    // fill index lists
    bool p1 = true;
    for(unsigned i = 0; i < this->algorithm.getGenomeCount(); i++) {

        // Check crossover indices to see if any match the current index
        for(unsigned j = 0; j < crossoverIndices.size(); j++)
            if(crossoverIndices[j] == i)
                p1 = !p1;

        // Add indices to the index lists
        if(i % 2) p1Indices.push_back(i);
        else p2Indices.push_back(i);
    }

    // If there are indices in group one
    if(p1Indices.size() != 0) {
        if(pop1.isLocal()) {    // Source population one is local
            if(this->isLocal()) this->copyGenomes(p1Indices, pop1);
            else pop1.exportGenomes(p1Indices, *this);
        } else {                // Source population is nonlocal
            if(this->isLocal()) this->importGenomes(p1Indices, pop1);
        }
    }

    // If there are indices in group two
    if(p2Indices.size() != 0) {
        if(pop2.isLocal()) {    // Source population two is local
            if(this->isLocal()) this->copyGenomes(p2Indices, pop2);
            else pop2.exportGenomes(p2Indices, *this);
        } else {                // Source population is nonlocal
            if(this->isLocal()) this->importGenomes(p2Indices, pop2);
        }
    }
}



// Get a copy of a specific genome
vector<genome> subPopulation::getGenomes(void) {

    // Check that this is a local subpopulation
    assertInitialised("Error, attempt to retrieve genomes from uninitialised subpopulation.");
    assertLocal("Error, attempt to retrieve genomes from nonlocal subpopulation.");

    // Return the genomes
    return this->genomes;
}



// Print out the subpopulation rankmap
void subPopulation::printRankMap(genomeTarget& target) {

    // If this isn't local, do nothing
    if(!this->isLocal()) {
        return;
    }

    cout << "Domain index: " << this->domainIndex << " on " << rankString() << "\n";
    for(unsigned i = 0; i < this->rankMap.size(); i++) {
        cout << rankMap[i].ptr->getPerfData(target).str() << "\n";
    }
}
