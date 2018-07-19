// Standard headers
#include "stdint.h"
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
using namespace std;


// Project headers
#include "bitVector.hpp"



// Basic empty constructor
bitVector::bitVector(void) {
    this->reset();
}


// Constructor which takes an argument for size
bitVector::bitVector(uint32_t l) {
    this->reset();
    this->init(l);
}


// Reset function, resets all internal state
void bitVector::reset(void) {
    this->length = 0;
    this->bitmaps.clear();
}


// Initialise the bit vector to given width
// Destroys present contents of the vector
void bitVector::init(uint32_t l) {

    // Reset everything
    this->length = l;
    this->bitmaps.clear();

    // Calculate number of bitmaps required to serve as base storage
    uint32_t bitmapCount = this->length / 64;
    if(this->length % 64) bitmapCount++;
    this->bitmaps.reserve(bitmapCount);
    for(unsigned i = 0; i < bitmapCount; i++) {
        this->bitmaps.push_back(0);
    }
}


// Returns bitmap index for a given bit number
uint32_t bitVector::bitmapIndex(uint32_t bitIndex) {
    return bitIndex / 64;
}


// Returns bitmap offset for a given bit number
uint32_t bitVector::bitmapOffset(uint32_t bitIndex) {
    return 63 - (bitIndex % 64);
}


// Returns mask for a given bitmap
uint64_t bitVector::bitmapMask(uint32_t bitmapIndex) {

    // Calculate maximum allowable index
    uint32_t maxIndex = bitmaps.size() - 1;
    uint64_t mask;

    // Return the apropriate mask
    if(bitmapIndex < maxIndex) {
        mask = ~((uint64_t)0);
    } else if(bitmapIndex == maxIndex) {
        uint32_t finalBitOffset = this->bitmapOffset(this->length - 1);
        mask = ~(((uint64_t)0x01 << finalBitOffset) - 1);
    } else {
        mask = 0;
    }

    // Return the mask
    return mask;
}


// Prints out bits in the bit vector
string bitVector::str(uint8_t format) {
    stringstream ss;

    // Switch by appropriate format
    switch(format) {

        // Hexadecimal
        case BITVECTOR_FORMAT_HEX:
            for(unsigned i = 0; i < this->bitmaps.size(); i++) {
                ss << hex << setfill('0') << setw(16) << (uint64_t)this->bitmaps[i];
            }
            break;

        // Binary
        case BITVECTOR_FORMAT_BIN:
            for(unsigned i = 0; i < this->length; i++) {
                if(this->getBit(i)) {
                    ss << "1";
                } else {
                    ss << "0";
                }
            }
            break;

        // Format not recognised
        default:
            return string("");
            break;
    }

    // Return the contents of the stringstream
    return ss.str();
}


// Get bitmap in the input bitmap vector
uint64_t bitVector::getBitmap(uint32_t bitmapIndex) {
    if(bitmapIndex >= this->bitmaps.size()) {
        cout << "Error, attempt to access bitmap beyond range.\n";
        exit(1);
    }
    return this->bitmaps[bitmapIndex];
}


// Set a bitmap in the input bitmap vector
void bitVector::setBitmap(uint32_t bitmapIndex, uint64_t bitmapValue) {
    if(bitmapIndex >= this->bitmaps.size()) {
        cout << "Error, attempt to access bitmap beyond range.\n";
        exit(1);
    }
    this->bitmaps[bitmapIndex] = bitmapValue & this->bitmapMask(bitmapIndex);
}


// Get specific bit
uint8_t bitVector::getBit(uint32_t bitIndex) {

    // Check that bit is within range of bit vector
    if(bitIndex >= this->length) {
        cout << "Error, attempt to read bit " << bitIndex << " from " << this->length << "-bit vector.\n";
        exit(1);
    }

    // Calculate bitmap index and offset
    uint32_t bitmapIndex = this->bitmapIndex(bitIndex);
    uint32_t bitmapOffset = this->bitmapOffset(bitIndex);

    // Retrieve apropriate bitmap
    uint64_t bitmap = this->bitmaps[bitmapIndex];
    if(bitmap & ((uint64_t)0x01 << bitmapOffset)) {
        return 1;
    } else {
        return 0;
    }
}


// Set specific bit
void bitVector::setBit(uint32_t bitIndex, uint8_t bitValue) {

    // Check that bit is within range of bit vector
    if(bitIndex >= this->length) {
        cout << "Error, attempt to write bit " << bitIndex << " in a " << this->length << "-bit vector.\n";
        exit(1);
    }

    // Calculate bitmap index and offset
    uint32_t bitmapIndex = this->bitmapIndex(bitIndex);
    uint32_t bitmapOffset = this->bitmapOffset(bitIndex);

    // Write the bit to the apropriate bitmap element
    if(bitValue == 0) {
        this->bitmaps[bitmapIndex] &= ~((uint64_t)0x01 << bitmapOffset);
    } else {
        this->bitmaps[bitmapIndex] |= (uint64_t)0x01 << bitmapOffset;
    }
}


// Append bit to end of vector
// zero = 0, nonzero = 1
void bitVector::appendBit(uint8_t bitValue) {

    // Calculate new length and number of required bitmaps
    this->length++;
    uint32_t requiredBitmaps = this->length / 64;
    if(this->length % 64) requiredBitmaps++;

    // Add a bitmap if neccessary
    if(this->bitmaps.size() < requiredBitmaps) {
        bitmaps.push_back(0);
    }

    // Write bit
    this->setBit(this->length - 1, bitValue);
}
