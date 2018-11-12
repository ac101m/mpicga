// Standard libs
#include <iostream>
#include <sstream>


// Project headers
#include "truthTable.hpp"


// Generates a multiplier
void generateMultiplier(string path, string widthStr) {

  // Get width of multiplier
  unsigned multiplierWidth;
  stringstream ss(widthStr);
  ss >> multiplierWidth;

  // Calculate input and output counts, create the truth table
  unsigned inputCount = multiplierWidth * 2;
  unsigned outputCount = multiplierWidth * 2;
  truthTable t(inputCount, outputCount);

  // Fill the truth table with the multiplier pattern
  unsigned patternCount = 0x01 << inputCount;
  for(unsigned i = 0; i < patternCount; i++) {
    unsigned a = i & ((0x01 << multiplierWidth) - 1);
    unsigned b = (i >> multiplierWidth) & ((0x01 << multiplierWidth) - 1);
    t.addPattern(make_pair(i, a * b));
  }

  // Output the pattern to the file
  t.writeToFile(path);
}


// Generate an adder
void generateAdder(string path, string widthStr, bool doCarry) {

  // Get width of adder
  unsigned adderWidth;
  stringstream ss(widthStr);
  ss >> adderWidth;

  // Determine input/output counts
  unsigned inputCount = adderWidth * 2;
  unsigned outputCount = adderWidth;
  if(doCarry) {
    inputCount++;
    outputCount++;
  }

  // Instantiate
  truthTable t(inputCount, outputCount);

  // Fill the truth table with the adder pattern
  unsigned patternCount = (0x01 << inputCount);
  for(unsigned i = 0; i < patternCount; i++) {
    unsigned a = i & ((0x01 << adderWidth) - 1);
    unsigned b = (i >> adderWidth) & ((0x01 << adderWidth) - 1);
    unsigned c = (i >> (adderWidth * 2)) & 0x01;
    t.addPattern(i, a + b + c);
  }

  // Write the pattern to the file
  t.writeToFile(path);
}


// Generates a pattern
int main(int argc, char** argv) {

  // Check that there are enough arguments
  if(argc < 3) {
    cout << "Usage: " << argv[0] << " [filename] [pattern] <pattern args>\n";
    cout << "\t available patterns: add, mul.\n";
    return 0;
  }

  // Process the arguments
  if(string(argv[2]) == "mul") {
    if(argc < 4) {
      cout << "Usage: " << argv[2] << " <input width>\n";
    } else {
      generateMultiplier(string(argv[1]), string(argv[3]));
    }

  } else if (string(argv[2]) == "add") {
    if(argc < 5) {
      cout << "Usage: " << argv[2] << " <input width> <carry=true>\n";
    } else {
      if(string(argv[4]) == "carry=true") {
        generateAdder(string(argv[1]), string(argv[3]), true);
      } else if(string(argv[4]) == "carry=false") {
        generateAdder(string(argv[1]), string(argv[3]), false);
      } else {
        cout << "Usage: " << argv[2] << " <input width> <carry=[true,false]>\n";
      }
    }

  } else {
    cout << "Unrecognised pattern: '" << argv[2] << "'\n";
  }

  // El fin
  return 0;
}
