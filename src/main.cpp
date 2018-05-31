// CPP
#include <iostream>
#include <random>
#include <functional>
using namespace std;


// Project related stuff
#include "config.hpp"
#include "utils.hpp"
#include "mpicga.hpp"
#include "bitVector.hpp"


// MPI
#include "mpi.h"


// Input width of multiplier pattern
//#define ADDER_INPUT_WIDTH 4
#define MULTIPLIER_INPUT_WIDTH 4



// Define the fitness function
uint32_t genomeFF(subPopulationPerf_t perf) {
    return perf.bestGenomeFitness;
}



// Define the fitness function
uint32_t subPopFF(genomePerf_t perf) {
    uint32_t effectiveActiveGenes = perf.activeGenes;
    if(perf.bitErrors) effectiveActiveGenes = 1024;
    return (perf.bitErrors << 6) + (effectiveActiveGenes << 3) + perf.genomeAge;
}



// Main routine
int main(int argc, char **argv) {

    // Initialise MPI
    MPI_Init(&argc, &argv);


    // Multiplier pattern
#ifdef MULTIPLIER_INPUT_WIDTH
    genomeTarget target(MULTIPLIER_INPUT_WIDTH * 2, MULTIPLIER_INPUT_WIDTH * 2);
    for(unsigned i = 0; i < 0x01 << (MULTIPLIER_INPUT_WIDTH * 2); i++) {

        // Calculate input and output words
        uint32_t mask = (0x01 << MULTIPLIER_INPUT_WIDTH) - 1;
        uint32_t a = i & mask;
        uint32_t b = (i >> MULTIPLIER_INPUT_WIDTH) & mask;

        // Add a pattern to the pattern map
        target.addPattern(i, a * b);
    }
#endif // MULTIPLIER


    // Adder pattern
#ifdef ADDER_INPUT_WIDTH
    genomeTarget target(ADDER_INPUT_WIDTH * 2, ADDER_INPUT_WIDTH + 1);
    for(unsigned i = 0; i < 0x01 << (ADDER_INPUT_WIDTH * 2); i++) {

        // Calculate input and output words
        uint32_t mask = (0x01 << ADDER_INPUT_WIDTH) - 1;
        uint32_t a = i & mask;
        uint32_t b = (i >> ADDER_INPUT_WIDTH) & mask;

        // Add a pattern to the pattern map
        target.addPattern(i, a + b);
    }
#endif // ADDER


    // Subpopulation distribution across ranks counts
    uint32_t subPopulationCount = 24;
    uint32_t totalGenerations = 32 * 32 * 480;
    uint32_t generationsPerSubPopulation = totalGenerations / subPopulationCount;
    uint32_t generationsPerCycle = 1024;
    uint32_t cycleCount = (totalGenerations / subPopulationCount) / generationsPerCycle;


    // zeroth rank, print out run information
    if(myRank() == 0) {
        cout << "[GENERATION CONFIG]\n";
        cout << "Total generations: " << totalGenerations << "\n";
        cout << "Generations per sub population: " << generationsPerSubPopulation << "\n";
        cout << "Generations per cycle: " << generationsPerCycle << "\n";
        cout << "Cycle count: " << cycleCount << "\n";
        cout << "\n[POPULATION PROCESS DISTRIBUTION]\n";
        cout << "Process count: " << rankCount() << "\n";
        cout << "Sub population count: " << subPopulationCount << "\n";
        cout << "Subpopulations per process: " << subPopulationCount / rankCount() << "\n";
        cout << "\n";
    }


    // Create a population and start timing
    population p(subPopulationCount, 8, 1024);


    // Population algorithm settings
    p.getAlgorithm().setGenerationsPerCycle(generationsPerCycle);
    p.getAlgorithm().setSeed(1);
    p.getAlgorithm().setCrossoverCount(3);
    p.getAlgorithm().setSelectCount(1);


    // Subpopulation algorithm settings
    p.getAlgorithm().getSubPopulationAlgorithm().setMutateCount(2);
    p.getAlgorithm().getSubPopulationAlgorithm().setMutateCount(1);
    p.getAlgorithm().getSubPopulationAlgorithm().setAllowableFunctions({GENE_FN_NAND});
    p.initialise(target, subPopFF);


    // Start ze timer
    double startTime = MPI_Wtime();


    // Iterate the population here
    p.iterate(target, subPopFF, cycleCount);


    // End ze timer
    double endTime = MPI_Wtime();


    // Quick barrier to stop execution duration overwriting stuff
    MPI_Barrier(MPI_COMM_WORLD);


    // Print ze timer
    if(myRank() == 0) {
        cout << "\nTotal execution time: " << endTime - startTime << "s\n";
    }

    // El fin
    MPI_Finalize();
    return 0;
}
