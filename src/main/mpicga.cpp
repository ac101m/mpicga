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
#include "optparse.hpp"



// Build option parser
OptionParser buildOptionParser(int argc, char **argv) {
  OptionParser options = OptionParser(
    argc, argv,
    "mpicga - a parallel genetic algorithm for generating cobinational logic circuits.");

  options.Add(Option("subpopcount", 'n', ARG_TYPE_INT,
                     "Set number of subpopulations for the algorithm to use.",
                     {DEFAULT_SUBPOP_COUNT}));

  options.Add(Option("subpopsize", 'S', ARG_TYPE_INT,
                     "Set size of subpopulations.",
                     {DEFAULT_SUBPOP_SIZE}));

  options.Add(Option("genomesize", 's', ARG_TYPE_INT,
                     "Set length of genomes.",
                     {DEFAULT_GENOME_SIZE}));

  options.Add(Option("totalgenerations", 'G', ARG_TYPE_INT,
                     "Set total number of generations for this run.",
                     {DEFAULT_TOTAL_GENERATIONS}));

  options.Add(Option("generationspercycle", 'g', ARG_TYPE_INT,
                     "Set number of generations per sub-population cycle.",
                     {DEFAULT_GENERATIONS_PER_CYCLE}));

  options.Add(Option("patternfile", 'p', ARG_TYPE_STRING,
                     "Path to file containing target pattern.",
                     {DEFAULT_PATTERN_PATH}));

  return options;
}


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
      while (0 == i) {
        sleep(5);
      }
    }
    #endif // DO_DEBUG_ATTACH

    // Build the option parser
    OptionParser options = buildOptionParser(argc, argv);

    // Load the pattern from file
    truthTable target(options.Get("patternfile"));

    // Subpopulation distribution across ranks counts
    int subPopCount = options.Get("subpopcount");
    int totalGenerations = options.Get("totalgenerations");
    int generationsPerCycle = options.Get("generationspercycle");
    int subPopSize = options.Get("subpopsize");
    int genomeSize = options.Get("genomesize");
    uint32_t generationsPerSubPopulation = totalGenerations / subPopCount;
    uint32_t cycleCount = (totalGenerations / subPopCount) / generationsPerCycle;

    // Initialise MPI
    MPI_Init(&argc, &argv);

    // zeroth rank, print out run information
    if(myRank() == 0) {
        cout << "\n[GENERATION CONFIG]\n";
        cout << "Total generations: " << totalGenerations << "\n";
        cout << "Generations per sub population: " << generationsPerSubPopulation << "\n";
        cout << "Generations per cycle: " << generationsPerCycle << "\n";
        cout << "Cycle count: " << cycleCount << "\n";
        cout << "\n[POPULATION LAYOUT]\n";
        cout << "Genome length: " << genomeSize << "\n";
        cout << "Subpopulation size: " << subPopSize << "\n";
        cout << "Total genomes: " << subPopCount * subPopSize << "\n";
        cout << "\n[PROCESS DISTRIBUTION]\n";
        cout << "Process count: " << rankCount() << "\n";
        cout << "Sub population count: " << subPopCount << "\n";
        cout << "Subpopulations per process: " << subPopCount / rankCount() << "\n\n";
    }

    // Create a population and start timing
    population p(subPopCount, subPopSize, genomeSize);

    // Population algorithm settings
    p.getAlgorithm().setGenerationsPerCycle(generationsPerCycle);
    p.getAlgorithm().setSeed(1);
    p.getAlgorithm().setCrossoverCount(3);
    p.getAlgorithm().setSelectCount(1);

    // Subpopulation algorithm settings
    p.getAlgorithm().getSubPopulationAlgorithm().setMutateCount(1);
    p.getAlgorithm().getSubPopulationAlgorithm().setAllowableFunctions({GENE_FN_NAND});
    p.initialise(target, subPopFF);

    // Iterate the population here
    double startTime = MPI_Wtime();
    p.iterate(target, subPopFF, cycleCount);
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
