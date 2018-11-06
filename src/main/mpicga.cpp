// Standard headers
#include <iostream>
#include <random>
#include <functional>
#include <unistd.h>
using namespace std;


// Project headers
#include "mpi.h"
#include "config.hpp"
#include "utils.hpp"
#include "mpicga.hpp"
#include "bitVector.hpp"



// Define the fitness function for subpopulations
uint32_t genomeFF(subPopulationPerf_t perf) {
    return perf.bestGenomeFitness;
}



// Define the fitness function for genomes
uint32_t subPopFF(genomePerf_t perf) {
    uint32_t effectiveActiveGenes = perf.activeGenes;
    if(perf.bitErrors) effectiveActiveGenes = 1024;
    return (perf.bitErrors << 6) + (effectiveActiveGenes << 3) + perf.genomeAge;
}



// Main routine
int main(int argc, char **argv) {


    // GDB attach point, for when shit gets squirly
    #ifdef DO_DEBUG_ATTACH
    if(DO_DEBUG_ATTACH) {
        int i = 0;
        char hostname[256];
        gethostname(hostname, sizeof(hostname));
        printf("PID %d on %s ready for attach\n", getpid(), hostname);
        fflush(stdout);
        while (0 == i)
            sleep(5);
    }
    #endif // DO_DEBUG_ATTACH


    // Initialise MPI
    MPI_Init(&argc, &argv);


    // Load the pattern from file
    truthTable target(DEFAULT_PATTERN_PATH);


    // Subpopulation distribution across ranks counts
    uint32_t subPopulationCount = 1920;
    uint32_t totalGenerations = 1024 * 1024 * 240;
    uint32_t generationsPerSubPopulation = totalGenerations / subPopulationCount;
    uint32_t generationsPerCycle = 256;
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
    population p(subPopulationCount, 4, 4096);


    // Population algorithm settings
    p.getAlgorithm().setGenerationsPerCycle(generationsPerCycle);
    p.getAlgorithm().setSeed(1);
    p.getAlgorithm().setCrossoverCount(3);
    p.getAlgorithm().setSelectCount(1);


    // Subpopulation algorithm settings
    p.getAlgorithm().getSubPopulationAlgorithm().setMutateCount(1);
    p.getAlgorithm().getSubPopulationAlgorithm().setAllowableFunctions({GENE_FN_AND, GENE_FN_OR, GENE_FN_XOR});
    p.initialise(target, subPopFF);


    // Start ze timer
    double startTime = MPI_Wtime();


    // Iterate the population here
    p.iterate(target, subPopFF, cycleCount);


    // End ze timer
    double endTime = MPI_Wtime();


    // Quick barrier to stop execution duration overwriting stuff
    MPI_Barrier(MPI_COMM_WORLD);


    // Print time difference
    if(myRank() == 0) {
        cout << "\nTotal execution time: " << endTime - startTime << "s\n";
    }

    // El fin
    MPI_Finalize();
    return 0;
}