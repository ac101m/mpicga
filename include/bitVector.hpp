#ifndef BITVECTOR_H
#define BITVECTOR_H


// Standard
#include "stdint.h"
#include <vector>
#include <string>


// Print formats
#define BITVECTOR_FORMAT_HEX 0
#define BITVECTOR_FORMAT_BIN 1


// Bit vector class
class bitVector {
  private:
    uint32_t length;
    std::vector<uint64_t> bitmaps;

    // Internal management functions
    uint32_t bitmapIndex(uint32_t bit);         // Indexing of bitmap vector
    uint32_t bitmapOffset(uint32_t bit);        // Bit indexing within bitmap

  public:
    // Constructors
    bitVector(void);
    bitVector(uint32_t l);

    // General utility
    void reset(void);
    void init(uint32_t l);
    std::string str(uint8_t format);

    // Low level access
    uint64_t getBitmap(uint32_t bitmapIndex);
    void setBitmap(uint32_t bitmapIndex, uint64_t bitmapValue);
    uint32_t getLength(void) {return this->length;}
    uint32_t getBitmapCount(void) {return this->bitmaps.size();}
    uint64_t bitmapMask(uint32_t bitmapIndex);

    // Normal access
    void appendBit(uint8_t bitValue);
    uint8_t getBit(uint32_t bitIndex);
    void setBit(uint32_t bitIndex, uint8_t bitValue);
};


#endif // BITVECTOR_H
