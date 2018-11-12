// C standard stuff
#include <exception>
#include <iostream>


// Project headers
#include "truthTable.hpp"


// Test suite configuration
#define CATCH_CONFIG_MAIN
#include "catch.hpp"



TEST_CASE("Truth table class test", "[truthTable]") {

  unsigned multiplierTestWidth = 4;

  SECTION("Multiplier pattern test") {

    unsigned multiplierWidth = multiplierTestWidth;
    unsigned inputCount = multiplierWidth * 2;
    unsigned outputCount = multiplierWidth * 2;
    unsigned patternCount = 0x01 << inputCount;

    vector<pair<uint32_t, uint32_t>> testPatterns;
    for(unsigned i = 0; i < patternCount; i++) {
      unsigned a = i & ((0x01 << multiplierWidth) - 1);
      unsigned b = (i >> multiplierWidth) & ((0x01 << multiplierWidth) - 1);
      testPatterns.push_back(make_pair(i, a * b));
    }

    truthTable t(inputCount, outputCount);


    for(unsigned i = 0; i < testPatterns.size(); i++) {
      t.addPattern(testPatterns[i]);
    }


    SECTION("Correct initialisation of input and output counts") {
      REQUIRE(t.getInputCount() == inputCount);
      REQUIRE(t.getOutputCount() == outputCount);
    }


    SECTION("Patterns are the same when read out as they were when written in") {
      unsigned errorCount = 0;
      for(unsigned i = 0; i < testPatterns.size(); i++)
        if(t.getPattern(i) != testPatterns[i])
          errorCount++;

      REQUIRE(errorCount == 0);
    }


    SECTION("Correct recall of saved table parameters") {
      REQUIRE(t.getInputCount() == inputCount);
      REQUIRE(t.getOutputCount() == outputCount);
    }


    SECTION("Correct recall of saved patterns") {
      t.writeToFile("test.pat");
      t = truthTable("test.pat");
      unsigned errorCount = 0;
      for(unsigned i = 0; i < testPatterns.size(); i++)
        if(t.getPattern(i) != testPatterns[i])
          errorCount++;

      REQUIRE(errorCount == 0);
    }
  }
}
