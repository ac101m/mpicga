#ifndef TRUTH_TABLE_HPP
#define TRUTH_TABLE_HPP


// Standard
#include <map>
#include <vector>
#include <string>
#include <fstream>


// Internal
#include "bitVector.hpp"


// Parse exceptions
struct truthTableParseException : public std::exception {
  std::string msg;
  truthTableParseException(std::string s) : msg(s) {}
  ~truthTableParseException() throw() {}
  const char* what() const throw() {return msg.c_str();}
};


// Truth table file pointer
class TTfp {
  private:

    std::ifstream fp;    // File pointer
    uint32_t line;
    uint32_t column;

  public:

    TTfp(std::string path);                 // Constructor, opens the file
    char current(void);                     // Get current character
    void advance(void);                     // Advance the file pointer by 1 byte
    void assertCurrent(char ch);            // Assert the current character
    void skipLine(void);                    // Skips to end of line
    void skip(std::string chars);           // Skip over chars of given type
    std::string get(std::string chars);     // Get a string made up of these chars
    uint32_t getNumber(uint32_t radix);     // Get number with an arbitrary radix

    // Routines for getting the things vector
    std::pair<uint32_t, uint32_t> getPattern(uint32_t radix);
    void getPatternList(std::vector<std::pair<uint32_t, uint32_t>>& patterns, uint32_t radix);

    // Generate line string to indicate where problems occur
    std::string lineString(void);

    // File format character strings
    static const std::string whiteSpaceChars;
    static const std::string numberChars;
    static const std::string nameChars;
};


// Class to represent input and output patterns for a genome
class truthTable {
  private:

    std::map<uint32_t, uint32_t> patternMap;     // Map of patterns, used for error checking
    std::vector<bitVector> inputs;               // Vectors containing input patterns
    std::vector<bitVector> outputs;              // Vectors containing output patterns

  public:     // Public interface

    truthTable(std::string path);                                // Constructor, from file
    truthTable(uint32_t inputCount, uint32_t outputCount);  // Constructor, empty
    void addPattern(uint32_t iPattern, uint32_t oPattern);  // Add a pattern to the target
    void addPattern(std::pair<uint32_t, uint32_t> pattern);      // Add a pattern to the target
    void assertValid(void);                                 // Errors if table is not valid

    // Get for input and output bit counts
    uint32_t getInputCount(void) {return this->inputs.size();}
    uint32_t getOutputCount(void) {return this->outputs.size();}

    // Gets for patterns and pattern count
    uint32_t getPatternCount(void) {return this->patternMap.size();}
    std::pair<uint32_t, uint32_t> getPattern(uint32_t index);

    // Gets and sets for various bitmap related stuff
    uint32_t getBitmapCount(void) {return this->inputs[0].getBitmapCount();}
    uint64_t getInputBitmap(uint32_t inputIndex, uint32_t bitmapIndex);
    uint64_t getOutputBitmap(uint32_t outputIndex, uint32_t bitmapIndex);
    uint64_t getBitmapMask(uint32_t bitmapIndex);

    // File writing routines
    void writeToFile(std::string path, uint32_t radix);
    void writeToFile(std::string path);
};


#endif // TRUTH_TABLE_HPP
