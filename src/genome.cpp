// Standard headers
#include <iostream>
#include <vector>
using namespace std;


// Project headers
#include "bitVector.hpp"
#include "utils.hpp"
#include "mpicga.hpp"



// Initialisation function
genome::genome(uint32_t geneCount, subPopulationAlgorithm& algorithm) {

    // Allocate vector for gene information
    this->genes.clear();
    for(unsigned i = 0; i < geneCount; i++) {
        this->genes.push_back(gene());
    }

    // Set up each gene
    for(unsigned i = 0; i < geneCount; i++) {

        // Randomly select a gene function
        this->genes[i].setFunction(algorithm.randomGeneFunction());

        // Randomly select input indices, don't do this for gene 0
        if(i) {
            this->genes[i].setAIndex(algorithm.randomGeneInputIndex(i));
            this->genes[i].setBIndex(algorithm.randomGeneInputIndex(i));
        }
    }

    // Clear genome performance data
    this->perfData.bitErrors = 0;
    this->perfData.activeGenes = 0;
    this->perfData.maxGateDelays = 0;
    this->perfData.genomeAge = 0;

    // This is the finished
    this->perfDataValid = false;
}



// Evaluates genome performance and applies returns performance info
void genome::updatePerfData(truthTable& target) {

    // Genome performance data
    this->perfData.bitErrors = 0;
    this->perfData.activeGenes = 0;
    this->perfData.maxGateDelays = 0;

    // Check that target has inputs and outputs
    target.assertValid();

    // Outer loop iterates over bitmaps
    for(unsigned i = 0; i < target.getBitmapCount(); i++) {

        // First loop iterates over genes, invalidating all of the output buffers
        for(unsigned j = 0; j < this->genes.size(); j++) {
            this->genes[j].outputBufferValid = false;
        }

        // Second loop reaplies inputs
        for(unsigned j = 0; j < target.getInputCount(); j++) {
            this->genes[j].overrideBuffer(target.getInputBitmap(j, i));
        }

        // Last loop calculates outputs for all output genes
        // Also calculates and sums bit errors.
        // k iterates over output genes
        // j iterates over output target pattern bitmaps
        uint32_t k = this->genes.size() - target.getOutputCount();
        for(unsigned j = 0; j < target.getOutputCount(); j++) {
            uint64_t buffer, difference;

            // Get the buffer for the first output gene
            buffer = this->genes[k].getOutputBuffer(this->genes);

            // Compare it to the target, calculate bit errors
            difference = buffer ^ target.getOutputBitmap(j, i);
            difference &= target.getBitmapMask(i);
            this->perfData.bitErrors += countBits(difference);

            // Increment k for next loop
            k++;
        }
    }

    // Iterate over genome, calculate number of active genes
    this->perfData.activeGenes = 0;
    for(unsigned i = target.getInputCount(); i < this->genes.size(); i++) {
        if(this->genes[i].isActive()) {
            this->perfData.activeGenes++;
        }
    }

    // Indicate that performance data is now valid
    this->perfDataValid = true;
}



// Get performance data
genomePerf_t genome::getPerfData(truthTable& target) {

    // Check if performance data is valid
    if(!this->perfDataValid) {
        this->updatePerfData(target);
    }

    // Return performance data
    return this->perfData;
}



// Mutate the genome
void genome::mutate(subPopulationAlgorithm& algorithm) {

    // Iterate
    for(unsigned i = 0; i < algorithm.getMutateCount(); i++) {

        // Select a gene at random to mutate
        int32_t selectedGeneIdx = algorithm.localRand(1, this->genes.size() - 1);

        // New mutate code
        if(this->genes[selectedGeneIdx].mutate(selectedGeneIdx, algorithm)) {
            this->perfDataValid = false;
        }
    }

    // Invalidate the performance data (assumes mutation generated bit errors)
    this->perfData.genomeAge = 0;
}



// Parse the genome from an array of gene network frames
void genome::parseGeneNetworkFrameArray(geneNetworkFrame_t *networkFrameArray) {

    // Iterate over all genes, intialising them with the network frames
    for(unsigned i = 0; i < this->genes.size(); i++) {
        this->genes[i] = gene(networkFrameArray[i]);
    }

    // Performance data is now invalid
    this->perfData.genomeAge = 0;
    this->perfDataValid = false;
}



// Copy gene data from another genome
void genome::copyFrom(genome& g) {

    // Iterate over genes
    for(unsigned i = 0; i < this->genes.size(); i++) {
        this->genes[i] = g.getGenes()[i];
    }

    // Reset perf-data
    this->perfData.genomeAge = 0;
    this->perfDataValid = false;
}
