#ifndef TRUTH_TABLE_HPP
#define TRUTH_TABLE_HPP


// Standard stuff
#include <map>
#include <vector>
#include <string>
#include <fstream>
using namespace std;


// My libs
#include "bitVector.hpp"


// Parse exceptions
struct truthTableParseException : public exception {
  string msg;
  truthTableParseException(string s) : msg(s) {}
  ~truthTableParseException() throw() {}
  const char* what() const throw() {return msg.c_str();}
};


// Truth table file pointer
class TTfp {
  private:

    ifstream fp;    // File pointer
    uint32_t line;
    uint32_t column;

  public:

    TTfp(string path);                      // Constructor, opens the file
    char current(void);                     // Get current character
    void advance(void);                     // Advance the file pointer by 1 byte
    void assertCurrent(char ch);            // Assert the current character
    void skipLine(void);                    // Skips to end of line
    void skip(string chars);                // Skip over chars of given type
    string get(string chars);               // Get a string made up of these chars
    uint32_t getNumber(uint32_t radix);     // Get number with an arbitrary radix

    // Routines for getting the things vector
    pair<uint32_t, uint32_t> getPattern(uint32_t radix);
    void getPatternList(vector<pair<uint32_t, uint32_t>>& patterns, uint32_t radix);

    // Generate line string to indicate where problems occur
    string lineString(void);

    // File format character strings
    static const string whiteSpaceChars;
    static const string numberChars;
    static const string nameChars;
};


// Class to represent input and output patterns for a genome
class truthTable {
  private:

    map<uint32_t, uint32_t> patternMap;     // Map of patterns, used for error checking
    vector<bitVector> inputs;               // Vectors containing input patterns
    vector<bitVector> outputs;              // Vectors containing output patterns

  public:     // Public interface

    truthTable(string path);                                // Constructor, from file
    truthTable(uint32_t inputCount, uint32_t outputCount);  // Constructor, empty
    void addPattern(uint32_t iPattern, uint32_t oPattern);  // Add a pattern to the target
    void addPattern(pair<uint32_t, uint32_t> pattern);      // Add a pattern to the target
    void assertValid(void);                                 // Errors if table is not valid

    // Get for input and output bit counts
    uint32_t getInputCount(void) {return this->inputs.size();}
    uint32_t getOutputCount(void) {return this->outputs.size();}

    // Gets for patterns and pattern count
    uint32_t getPatternCount(void) {return this->patternMap.size();}
    pair<uint32_t, uint32_t> getPattern(uint32_t index);

    // Gets and sets for various bitmap related stuff
    uint32_t getBitmapCount(void) {return this->inputs[0].getBitmapCount();}
    uint64_t getInputBitmap(uint32_t inputIndex, uint32_t bitmapIndex);
    uint64_t getOutputBitmap(uint32_t outputIndex, uint32_t bitmapIndex);
    uint64_t getBitmapMask(uint32_t bitmapIndex);

    // File writing routines
    void writeToFile(string path, uint32_t radix);
    void writeToFile(string path);
};


#endif // TRUTH_TABLE_HPP
