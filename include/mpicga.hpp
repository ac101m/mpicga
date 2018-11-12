#ifndef MPIGA_H
#define MPIGA_H


// Standard stuff
#include "stdint.h"
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <random>
#include <functional>
using namespace std;


// Project headers
#include "mpi.h"
#include "truthTable.hpp"


// Class pre-declarations
class populationAlgorithm;
class subPopulationAlgorithm;



//========[GENE]=================================================================================//

// Gene values
#define GENE_FN_NOP '-'
#define GENE_FN_NOT '~'
#define GENE_FN_AND '&'
#define GENE_FN_NAND '@'
#define GENE_FN_OR '|'
#define GENE_FN_NOR ':'
#define GENE_FN_XOR '^'
#define GENE_FN_XNOR 'v'



// Minimal gene datastructure for transmission over the network
typedef struct {
    char geneFunction;
    uint16_t aInputIndex;
    uint16_t bInputIndex;
} geneNetworkFrame_t;



// Gene class, represents a gene within a genome
class gene {
    public:

        bool outputBufferValid;         // Flag to indicate wether or not output value buffer is valid
        char geneFunction;              // Character representative of gene logic function
        uint16_t aInputIndex;           // Input index A
        uint16_t bInputIndex;           // Input index B
        uint64_t outputBuffer;          // Output buffer bits

    public:

        // Constructors
        gene(void);                 // Blank
        gene(geneNetworkFrame_t);   // From network frame
        uint64_t getOutputBuffer(vector<gene>& genes);          // Get gene output buffer recursively
        uint64_t computeBufferValue(uint64_t a, uint64_t b);    // Calculate gene output with given inputs

        // Gets and sets for gene functions
        char getGeneFunction(void) {return this->geneFunction;}
        void setGeneFunction(char gf) {this->geneFunction = gf;}
        bool isActive(void) {return this->outputBufferValid;}

        // Gets and sets for output buffer
        void clearBuffer(void) {this->outputBufferValid = false;}
        void overrideBuffer(uint64_t iv);

        // General operations
        void setFunction(char fn) {this->geneFunction = fn;}
        void setAIndex(uint32_t a) {this->aInputIndex = a;}
        void setBIndex(uint32_t b) {this->bInputIndex = b;}
        bool mutate(uint32_t minIndex, uint32_t maxIndex, vector<char> fnPool);
        bool mutate(uint32_t selectedIndex, subPopulationAlgorithm& algorithm);

        // Generates and returns a gene template
        geneNetworkFrame_t getNetworkFrame(void);
};



//========[GENOME]===============================================================================//

// Struct to contain genome performance data
typedef struct {

    uint32_t genomeAge;
    uint32_t bitErrors;
    uint32_t activeGenes;
    uint32_t maxGateDelays;

    // Produces printout string for genome performance struct
    string str(void) {
        stringstream ss;
        ss << this->bitErrors << " \t";
        ss << this->activeGenes << " \t";
        ss << this->maxGateDelays << " \t";
        ss << this->genomeAge;
        return ss.str();
    }
} genomePerf_t;



// Genome class, represents an individual
class genome {
    private:

        // Genome state information
        vector<gene> genes;                 // Gene array

        // Genome performance data relative to input pattern used during evaluation
        genomePerf_t perfData;
        bool perfDataValid;

    private:

        // Update performance data
        void updatePerfData(truthTable& target);

    public:

        // Constructor
        genome(uint32_t geneCount, subPopulationAlgorithm& algorithm);

        // Gets for genome data
        uint32_t getGeneCount(void) {return this->genes.size();}
        vector<gene> getGenes(void) {return this->genes;}
        genomePerf_t getPerfData(truthTable& target);

        // Operators
        void mutate(subPopulationAlgorithm& behaviour);
        void incrementAge(void) {this->perfData.genomeAge++;}

        // Parse genome from an array of genome network frames
        void parseGeneNetworkFrameArray(geneNetworkFrame_t *networkFrameArray);

        // Copy gene data from one genome to this one
        void copyFrom(genome& g);
};



//========[GENOME TRANSMIT BUFFER]===============================================================//

// Class contains a buffer
class genomeTransmissionBuffer {
    private:

        // Raw data storage buffer
        geneNetworkFrame_t *buffer;
        uint32_t currentGenes;
        uint32_t maxGenes;

    private:

        // Append a gene to the buffer
        void append(geneNetworkFrame_t g);

        // Append a gene to the buffer
        void append(gene g);

    public:

        // Constructor
        genomeTransmissionBuffer(uint32_t bufferLength);

        // Get the raw data
        geneNetworkFrame_t *getData(void);

        // Append a genome to the buffer
        void append(genome g);

        // Transmit the buffer
        void transmit(int32_t destination, int32_t tag);

        // Recieve an incoming buffer
        void receive(int32_t source, int32_t tag);

        // Destructor
        ~genomeTransmissionBuffer(void);
};



//========[SUB POPULATION ALGORITHM]=============================================================//

// Structure to contain population evolution specifications
// Class dictates the behaviour of a population during evolution
class subPopulationAlgorithm {
    private:

        // Subpopulation geometry
        uint32_t genomeCount;
        uint32_t genomeLength;

        // sub population algorithm stuff
        uint32_t selectCount;
        uint32_t lowSelectRange;
        uint32_t highSelectRange;

        // genome algorithm spec stuff
        uint32_t mutateCount;
        uint32_t minFeedForward;
        uint32_t maxFeedForward;
        vector<char> allowableFunctions;

        // Local random number generator
        mt19937 localRandEngine;

        // Current tag for transmit operations
        uint32_t txTagCounter;

    public:

        // constructors default and non-default
        subPopulationAlgorithm(void);
        subPopulationAlgorithm(uint32_t genomeCount, uint32_t genomeLength);

        // Gets for population geometry
        // No sets as these are only initialisable at initialisation of class
        uint32_t getGenomeCount(void) {return this->genomeCount;}
        uint32_t getGenomeLength(void) {return this->genomeLength;}

        // Get for high select fraction
        double getHighSelectRange(void) {return this->highSelectRange;}
        double getLowSelectRange(void) {return this->lowSelectRange;}

        // Get and set for selection count
        uint32_t getSelectCount(void) {return this->selectCount;}
        void setSelectCount(uint32_t sc) {this->selectCount = sc;}

        // Get and set for min gate delays (max feed forward)
        uint32_t getMaxFeedForward(void) {return this->maxFeedForward;}
        void setMaxFeedForward(uint32_t ff) {this->maxFeedForward = ff;}
        void setMinGateDelays(uint32_t gd);

        // Get and set for max gate delays (min feed forward)
        uint32_t getMinFeedForward(void) {return this->minFeedForward;}
        void setMinFeedForward(uint32_t ff) {this->minFeedForward = ff;}
        void setMaxGateDelays(uint32_t gd);

        // Get and set for mutation count
        uint32_t getMutateCount(void) {return this->mutateCount;}
        void setMutateCount(uint32_t mc) {this->mutateCount = mc;}

        // Get and set for allowable functions
        vector<char> getAllowableFunctions(void) {return this->allowableFunctions;}
        void setAllowableFunctions(vector<char> af) {this->allowableFunctions = af;}

        // Local random number generator
        int32_t localRand(int32_t minimum, int32_t maximum);
        void setSeed(uint32_t seed) {this->localRandEngine.seed(seed);}

        // Select high and low genome
        int32_t randomHighGenome(void);             // Random index for high genome
        int32_t randomLowGenome(void);              // Random index for low genome
        uint32_t randomGeneInputIndex(int32_t i);   // Random gene input index for gene i
        char randomGeneFunction(void);              // Random allowable gene function
};



//========[SUBPOPULATION]========================================================================//

// Subpopulation performance data struct
typedef struct {
    uint32_t bestGenomeFitness;
} subPopulationPerf_t;



// Datastructure for entries in genome rankmap
typedef struct {
    genome *ptr;
    uint32_t index;
    uint32_t fitness;
} genomeFitnessMapping_t;



// Class to represent subpopulation of genomes
class subPopulation {
    private:

        // Initialisation and rank information
        int32_t commWorldAddress;
        uint32_t domainIndex;
        bool initialised;
        bool local;

        // Genome mutation information
        subPopulationAlgorithm algorithm;

        // Population state data
        vector<genome> genomes;                     // Raw genome data
        vector<genomeFitnessMapping_t> rankMap;     // Genome rank map

    private:

        // Stuff for sorting the rankmap
        void swapRankMap(uint32_t i1, uint32_t i2);
        int64_t rankMapSortKey(uint32_t i);
        int32_t partitionRankMap(int32_t low, int32_t high);
        void quickSortRankMap(int32_t low, int32_t high);
        void sortRankMap(void);

        // Assert checks
        void assertInitialised(string msg);
        void assertLocal(string msg);

        // Internal copy transmit and recieve routines
        void parseGenomeBuffer(genomeTransmissionBuffer& buffer, vector<uint32_t>& genomeIndices);
        void copyGenomes(vector<uint32_t>& genomeIndices, subPopulation& source);
        void importGenomes(vector<uint32_t>& genomeIndices, subPopulation& source);
        void exportGenomes(vector<uint32_t>& genomeIndices, subPopulation& target);

    public:

        // Domain decomposition related
        uint32_t getDomainIndex(void) {return this->domainIndex;}
        int32_t getProcessRank(void) {return this->commWorldAddress;}

        // Constructor
        subPopulation(subPopulationAlgorithm algorithm);
        subPopulation(uint32_t populationSize, uint32_t genomeSize);

        // Get reference to algorithm datastructure
        subPopulationAlgorithm& getAlgorithm(void) {return this->algorithm;}

        // Initialise routine
        void initialise(truthTable& target, uint32_t(*ff)(genomePerf_t));
        void initialise(truthTable& target, uint32_t(*ff)(genomePerf_t), int32_t commWorldAddress);

        // Returns true if caller is the local process
        bool isLocal(void);

        // Iterate the population using specific mutation specs
        void iterate(truthTable& target, uint32_t(*ff)(genomePerf_t));
        void iterate(truthTable& target, uint32_t(*ff)(genomePerf_t), uint32_t n);

        // Get subpopulation performance data
        subPopulationPerf_t getPerfData(void);

        // Subpopulation crossover operator
        void crossover(subPopulation& pop1, subPopulation& pop2, vector<uint32_t> crossoverIndices);

        // Get a specific genome
        vector<genome> getGenomes(void);

        // Print out the rankmap
        void printRankMap(truthTable& target);

        // Update the rankmap
        void updateRankMap(truthTable& target, uint32_t(*ff)(genomePerf_t));
};



//========[POPULATION ALGORITHM]=================================================================//

// Class contains the algorithm specification for an entire population
// Specifies the behaviour of the population and all resident subpopulations
class populationAlgorithm {
    private:

        // Local random number generation
        mt19937 localRandEngine;

        // Crossover variables
        uint32_t selectCount;
        uint32_t crossoverCount;
        uint32_t lowSelectRange;
        uint32_t highSelectRange;

        // Subpopulation generations per synchronisation cycle
        uint32_t generationsPerCycle;

        // Subpopulation count algorithm and fitness function
        subPopulationAlgorithm subPopAlgorithm;
        uint32_t subPopulationCount;
        function<uint32_t(subPopulationPerf_t)> subPopulationFitnessFunction;

        // Processing modifiers
        int threadCount;

    public:

        // Constructors
        populationAlgorithm(void);
        populationAlgorithm(uint32_t subPopulationCount, uint32_t genomeCount, uint32_t genomeLength);

        // Get subpopulation algorithm
        subPopulationAlgorithm& getSubPopulationAlgorithm(void) {return this->subPopAlgorithm;}

        // Get subpopulation count
        uint32_t getSubPopulationCount(void) {return this->subPopulationCount;}

        // Get and set for generations per cycle
        uint32_t getGenerationsPerCycle(void) {return this->generationsPerCycle;}
        void setGenerationsPerCycle(uint32_t gpc) {this->generationsPerCycle = gpc;}

        // Local random number generator
        int32_t localRand(int32_t minimum, int32_t maximum);
        void setSeed(uint32_t seed) {this->localRandEngine.seed(seed);}

        // Get and set for crossover count
        uint32_t getCrossoverCount(void) {return this->crossoverCount;}
        void setCrossoverCount(uint32_t cc) {this->crossoverCount = cc;}

        // Get and set for crossover count
        uint32_t getSelectCount(void) {return this->selectCount;}
        void setSelectCount(uint32_t sc) {this->selectCount = sc;}

        // Get and set for thread count
        int getThreadCount(void) {return this->threadCount;}
        void setThreadCount(int tc) {this->threadCount = tc;}

        // Select subpopulations
        int32_t randomLowSubPopulation(void);
        int32_t randomHighSubPopulation(void);

        // Get a random list of crossover indices
        vector<uint32_t> randomCrossoverIndices(void);
};



//========[POPULATION]===========================================================================//

// Domain decomposition function
// Used to distribute populations across ranks
int32_t domainDecomposition(uint32_t indexWithinDomain);



// Datastructure for entries in subPopulation rankmap
typedef struct {
    subPopulation *ptr;
    uint32_t index;
    uint32_t fitness;
} subPopulationFitnessMapping_t;



// Class to represent whole domain
// This is a parallel class, data resides on multiple nodes
class population {
    private:

        // Population algorithm
        populationAlgorithm algorithm;

        // Is this initialised?
        bool initialised;

        // Vector of subpopulations
        vector<subPopulation> subPopulations;
        vector<subPopulationFitnessMapping_t> rankMap;
        vector<uint32_t> rankSubPopulationCounts;

        // Communication tag counter used to uniquely identify
        uint32_t commTagCounter;

    private:

        // Error and end if this is not initialised
        void assertInitialised(string msg);

        // Get local population
        vector<uint32_t> getLocalSubPopulationIndices(void);

        // Iterate the population n generations
        void iterateSubPopulations(truthTable& target, uint32_t(*ff)(genomePerf_t), uint32_t n);

        // Rank map sorting
        void swapRankMap(uint32_t i1, uint32_t i2);
        int64_t rankMapSortKey(uint32_t i);
        int32_t partitionRankMap(int32_t low, int32_t high);
        void quickSortRankMap(int32_t low, int32_t high);
        void sortRankMap(void);

        // Rankmap synchonisation and helper routines
        unsigned *getRankMapTxBuffer(void);
        void parseRankMapRxBuffer(unsigned *rxBuffer);
        void synchroniseRankMap(void);
        void updateRankMap(void);

        // Get local subpopulation counts
        uint32_t getLocalSubPopulationCount(void);              // From my rank
        uint32_t getSubPopulationCount(uint32_t rankAddress);   // From specific rank

        // Perform single crossover event
        void doSubPopulationCrossover(truthTable& target, uint32_t(*ff)(genomePerf_t));

    public:

        // Default constructor
        population(uint32_t subPopulationCount, uint32_t genomeCount, uint32_t genomeLength);

        // Get population algorithm
        populationAlgorithm& getAlgorithm(void) {return this->algorithm;}

        // Initialise
        void initialise(truthTable& target, uint32_t(*ff)(genomePerf_t));

        // Iterate the population using specific mutation specs
        void iterate(truthTable& target, uint32_t(*ff)(genomePerf_t));
        void iterate(truthTable& target, uint32_t(*ff)(genomePerf_t), uint32_t n);

        // Print the subpopulation rankmap
        void printRankMap(void);
};



#endif // MPIGA_H
