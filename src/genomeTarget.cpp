// Standard headers
#include <iostream>


// Project headers
#include "mpicga.hpp"



// Constructor
genomeTarget::genomeTarget(uint32_t inputCount, uint32_t outputCount) {

    // Check that input and count is valid
    if(inputCount == 0) {
        cout << "Error, input count must be nonzero.\n";
        exit(1);
    }

    // Check that output count is valid
    if(outputCount == 0) {
        cout << "Error, output count must be nonzero.\n";
        exit(1);
    }

    // Clear the input and output vectors
    this->inputs.clear();
    this->outputs.clear();

    // Reserve space for inputs
    this->inputs.reserve(inputCount);
    for(unsigned i = 0; i < inputCount; i++) {
        this->inputs.push_back(bitVector(0));
    }

    // Reserve space for outputs
    this->outputs.reserve(outputCount);
    for(unsigned i = 0; i < outputCount; i++) {
        this->outputs.push_back(bitVector(0));
    }
}


// Returns true of the pattern is a valid optimisation target
void genomeTarget::assertValid(void) {

    // First, get length of first input vector
    uint32_t bitPatternCount = this->inputs[0].getLength();

    // Compare it to other input vectors
    for(unsigned i = 1; i < this->inputs.size(); i++) {
        if(this->inputs[i].getLength() != bitPatternCount) {
            cout << "Error, mismatched input vector length in target pattern.\n";
            exit(1);
        }
    }

    // Compare it to output vectors
    for(unsigned i = 0; i < this->outputs.size(); i++) {
        if(this->outputs[i].getLength() != bitPatternCount) {
            cout << "Error, mismatched output vector length in target pattern.\n";
            exit(1);
        }
    }

    // Check that pattern vector contains patterns
    if(bitPatternCount == 0) {
        cout << "Error, target pattern is empty.\n";
        exit(1);
    }

    // Check that
    if(!this->inputs.size()) {
        cout << "Error, target pattern contains no input vectors.\n";
        exit(1);
    }

    // Check that
    if(!this->outputs.size()) {
        cout << "Error, target pattern contains no output vectors.\n";
        exit(1);
    }
}


// Add a pattern to the target
// LSB = input 0, and so forth.
void genomeTarget::addPattern(uint32_t iPattern, uint32_t oPattern) {

    // Mask the input and output patterns
    uint32_t iPatternMasked = iPattern & ((((uint32_t)0x01) << this->inputs.size()) - 1);
    uint32_t oPatternMasked = oPattern & ((((uint32_t)0x01) << this->inputs.size()) - 1);

    // Check whether input pattern has already been specified
    if(this->patternMap.find(iPatternMasked) != this->patternMap.end()) {
        if(this->patternMap[iPatternMasked] != oPatternMasked) {
            cout << "Error, conflicting target pattern specified, exiting.\n";
            exit(1);
        } else {
            cout << "Warning, duplicate pattern definition ignored\n";
            return;
        }
    }

    // Add pattern to input vectors
    for(unsigned i = 0; i < this->inputs.size(); i++) {
        if(iPatternMasked & (0x01 << i)) {
            this->inputs[i].appendBit(1);
        } else {
            this->inputs[i].appendBit(0);
        }
    }

    // Add to output vectors
    for(unsigned i = 0; i < this->outputs.size(); i++) {
        if(oPatternMasked & (0x01 << i)) {
            this->outputs[i].appendBit(1);
        } else {
            this->outputs[i].appendBit(0);
        }
    }

    // Add new pattern to pattern map
    this->patternMap.insert(make_pair(iPatternMasked, oPatternMasked));
}


// Gets the input bitmap associated with
uint64_t genomeTarget::getInputBitmap(uint32_t inputIndex, uint32_t bitmapIndex) {
    return this->inputs[inputIndex].getBitmap(bitmapIndex);
}


// Gets the output bitmap associated with
uint64_t genomeTarget::getOutputBitmap(uint32_t outputIndex, uint32_t bitmapIndex) {
    return this->outputs[outputIndex].getBitmap(bitmapIndex);
}


//
uint64_t genomeTarget::getBitmapMask(uint32_t bitmapIndex) {
    return this->inputs[0].bitmapMask(bitmapIndex);
}
