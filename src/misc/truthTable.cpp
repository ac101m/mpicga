// Standard headers
#include <iostream>
#include <sstream>
using namespace std;


// Project headers
#include "truthTable.hpp"


//========[FILE PARSING CLASS]===================================================================//


// Static pars
const string TTfp::whiteSpaceChars = " \t\n";
const string TTfp::numberChars = "0123456789abcdefABCDEF";
const string TTfp::nameChars = "ABCDEFGHIJKLMNOPGRSTUVWXYZabcdefghijklmnopqrstuvwxyz_0123456789";


// Constructor, opens the file
TTfp::TTfp(string path) {

    // Open the file for reading
    this->fp.open(path);

    // Check that the file is indeed open
    if(!fp.is_open()) {
        throw(runtime_error("Could not open file '" + path + "' for reading."));
    }

    // Reset line and column info
    this->line = 1;
    this->column = 0;
}


// Return the currently pointed to character
char TTfp::current(void) {
    return this->fp.peek();
}


// Advance the pointer
void TTfp::advance(void) {
    if(this->fp.get() == '\n') {
        this->column = 0;
        this->line++;
    } else {
        this->column++;
    }
}


// Assert the value of the current character
void TTfp::assertCurrent(char ch) {
    if(this->current() != ch) {
        stringstream ss;
        ss << "Unexpected '" << this->current();
        ss << "', expected '" << ch << "'.";
        throw(truthTableParseException(ss.str()));
    }
}


// Skip whitespace
void TTfp::skipLine(void) {
    do {
        if(this->current() == '\n') {
            this->advance();
            return;
        } else {
            this->advance();
        }
    } while(this->current() != EOF);
}


// Skip over chars found in a certain string
void TTfp::skip(string chars) {
    do {
        if(chars.find(this->current()) != string::npos) {
            this->advance();
        } else {
            return;
        }
    } while(this->current() != EOF);
}


// Get a sequence of characters from a certain subset of characters
string TTfp::get(string chars) {
    string s = "";
    while(chars.find(this->current()) != string::npos) {
        s += this->current();
        this->advance();
    }
    return s;
}


// Get bit pattern represented by number with arbitrary radix
uint32_t TTfp::getNumber(uint32_t radix) {

    // Make sure radix is apropriate
    if(radix < 2 || radix > 16) {
        throw(truthTableParseException("Supported radix values: 2 - 16."));
    }

    // Skip whitespace
    this->skip(whiteSpaceChars);

    // Get number string
    string text = this->get(TTfp::numberChars);

    // build a list of characters accepted at the specified radix
    string radixChars = "";
    for(unsigned i = 0; i < radix; i++) {
        radixChars += numberChars[i];
    }

    // Iterate
    uint32_t value = 0;
    for(unsigned i = 0; i < text.length(); i++) {
        if(radixChars.find(text[i]) == string::npos) {
            string ch = ""; ch += text[i];
            throw(truthTableParseException("'" + ch + "' outside radix bounds."));
        }
        switch(text[i]) {
            case '0': value *= radix; break;
            case '1': value *= radix; value += 1; break;
            case '2': value *= radix; value += 2; break;
            case '3': value *= radix; value += 3; break;
            case '4': value *= radix; value += 4; break;
            case '5': value *= radix; value += 5; break;
            case '6': value *= radix; value += 6; break;
            case '7': value *= radix; value += 7; break;
            case '8': value *= radix; value += 8; break;
            case '9': value *= radix; value += 9; break;
            case 'a': case 'A': value *= radix; value += 10; break;
            case 'b': case 'B': value *= radix; value += 11; break;
            case 'c': case 'C': value *= radix; value += 12; break;
            case 'd': case 'D': value *= radix; value += 13; break;
            case 'e': case 'E': value *= radix; value += 14; break;
            case 'f': case 'F': value *= radix; value += 15; break;
            default: return value;
        }
    }

    // Return the computed output
    return value;
}


// Parse individual pattern
pair<uint32_t, uint32_t> TTfp::getPattern(uint32_t radix) {
    uint32_t inputBitPattern, outputBitPattern;

    // Skip leading whitespace
    this->skip(whiteSpaceChars);

    // Start by reading a number at the given radix
    if(numberChars.find(this->current()) != string::npos) {
        inputBitPattern = this->getNumber(radix);
    } else {
        stringstream ss;
        ss << "Unexpected '" << this->current();
        ss << "', expected bit pattern specification.";
        throw(truthTableParseException(ss.str()));
    }

    // Skip more whitespace, look for divider
    this->skip(whiteSpaceChars);
    this->assertCurrent(':');
    this->advance();
    this->skip(whiteSpaceChars);

    // Start by reading a number at the given radix
    if(numberChars.find(this->current()) != string::npos) {
        outputBitPattern = this->getNumber(radix);
    } else {
        stringstream ss;
        ss << "Unexpected '" << this->current();
        ss << "', expected bit pattern specification.";
        throw(truthTableParseException(ss.str()));
    }

    // Add parsed pattern to the vector
    return pair<uint32_t, uint32_t>(inputBitPattern, outputBitPattern);
}


// Parse list of patterns
void TTfp::getPatternList(vector<pair<uint32_t, uint32_t>>& patternList, uint32_t radix) {

    do {
        // Parse a pattern
        patternList.push_back(this->getPattern(radix));

        // Skip trailing whitespace
        this->skip(whiteSpaceChars);

        // Check for either comma or semicolon
        if(this->current() == ';') {
            this->advance();
            return;
        } else if(this->current() == ',') {
            this->advance();
        } else {
            stringstream ss;
            ss << "Unexpected '" << this->current();
            ss << "', expected ';' or ','.";
            throw(truthTableParseException(ss.str()));
        }

    } while(this->current() != EOF);
}


// Generate line string for debugging
string TTfp::lineString(void) {
    stringstream ss;
    ss << "[Line " << this->line << ", col " << this->column << "]";
    return ss.str();
}



//========[PUBLIC INTERFACE METHODS]=============================================================//


// Constructor from file
truthTable::truthTable(string path) {

    // Parsing variables
    uint32_t radix = 0;
    uint32_t inputCount = 0;
    uint32_t outputCount = 0;
    vector<pair<uint32_t, uint32_t>> patterns;

    // Create a specialised parser/function pointer and skip leading whitespace
    TTfp fp(path);
    fp.skip(TTfp::whiteSpaceChars);

    // Parse the file and populate the table
    do {
        // Skip comments
        if(fp.current() == '#') {
            fp.skipLine();

        // Parse identifiers
        } else if(TTfp::nameChars.find(fp.current()) != string::npos) {
            string ident = fp.get(TTfp::nameChars);
            if(ident == "inputCount" || ident == "iCount") {
                if(inputCount) {
                    throw(truthTableParseException("Input count already specified."));
                } else {
                    inputCount = fp.getNumber(10);
                    fp.skip(TTfp::whiteSpaceChars);
                    fp.assertCurrent(';');
                    fp.advance();
                }

            } else if(ident == "outputCount" || ident == "oCount") {
                if(outputCount) {
                    throw(truthTableParseException("Output count already specified."));
                } else {
                    outputCount = fp.getNumber(10);
                    fp.skip(TTfp::whiteSpaceChars);
                    fp.assertCurrent(';');
                    fp.advance();
                }

            } else if(ident == "radix") {
                radix = fp.getNumber(10);
                fp.skip(TTfp::whiteSpaceChars);
                fp.assertCurrent(';');
                fp.advance();

            } else if(ident == "pattern") {
                if(radix) {
                    fp.getPatternList(patterns, radix);
                } else throw(truthTableParseException("Radix not specified."));

            } else {
                throw(truthTableParseException("Identifier '" + ident + "' not recognised."));
            }

        // Handle unexpected characters
        } else {
            stringstream ss;
            ss << fp.lineString() << " Unexpected '" << fp.current() << "'.";
            throw(truthTableParseException(ss.str()));
        }

        // Skip trailing whitespace
        fp.skip(TTfp::whiteSpaceChars);

    } while(fp.current() != EOF);

    // Check that input and output counts have been specified
    if(!inputCount) throw(truthTableParseException("File '" + path + "', input count not specified."));
    if(!outputCount) throw(truthTableParseException("File '" + path + "', output count not specified."));
    if(!patterns.size()) throw(truthTableParseException("File '" + path + "', no patterns specified."));

    // Initialise the table data structure and add all discovered patterns
    *this = truthTable(inputCount, outputCount);
    for(unsigned i = 0; i < patterns.size(); i++) {
        addPattern(patterns[i].first, patterns[i].second);
    }
}


// Constructor
truthTable::truthTable(uint32_t inputCount, uint32_t outputCount) {

    // Check that input and count is valid
    if(inputCount == 0) {
        throw(invalid_argument("Input count must be nonzero."));
    }

    // Check that output count is valid
    if(outputCount == 0) {
        throw(invalid_argument("Output count must be nonzero."));
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
void truthTable::assertValid(void) {

    // First, get length of first input vector
    uint32_t bitPatternCount = this->inputs[0].getLength();

    // Compare it to other input vectors
    for(unsigned i = 1; i < this->inputs.size(); i++) {
        if(this->inputs[i].getLength() != bitPatternCount) {
            throw(logic_error("Truth table consistency fail, input vector length mismatch."));
        }
    }

    // Compare it to output vectors
    for(unsigned i = 0; i < this->outputs.size(); i++) {
        if(this->outputs[i].getLength() != bitPatternCount) {
            throw(logic_error("Truth table consistency fail, output vector length mismatch."));
        }
    }

    // Check that
    if(!this->inputs.size()) {
        throw(logic_error("Truth table consistency fail, table contains no input vectors."));
    }

    // Check that
    if(!this->outputs.size()) {
        throw(logic_error("Truth table consistency fail, table contains no output vectors."));
    }

    // Check that pattern vector contains patterns
    if(bitPatternCount == 0) {
        throw(logic_error("Truth table consistency fail, table is empty."));
    }
}


// Add a pattern to the target
// LSB = input 0, and so forth.
void truthTable::addPattern(uint32_t iPattern, uint32_t oPattern) {

    // Mask the input and output patterns
    uint32_t iPatternMasked = iPattern & ((((uint32_t)0x01) << this->inputs.size()) - 1);
    uint32_t oPatternMasked = oPattern & ((((uint32_t)0x01) << this->outputs.size()) - 1);

    // Check whether input pattern has already been specified
    if(this->patternMap.find(iPatternMasked) != this->patternMap.end()) {
        if(this->patternMap[iPatternMasked] != oPatternMasked) {
            throw(logic_error("Truth table logic fail, conflicting pattern submitted."));
        } else {
            cout << "Warning, duplicate pattern [";
            cout << iPattern << ":" << oPattern << "], definition ignored\n";
            return;
        }
    }

    // Add pattern to input vectors
    for(unsigned i = 0; i < this->getInputCount(); i++) {
        if(iPatternMasked & (0x01 << i)) {
            this->inputs[i].appendBit(1);
        } else {
            this->inputs[i].appendBit(0);
        }
    }

    // Add to output vectors
    for(unsigned i = 0; i < this->getOutputCount(); i++) {
        if(oPatternMasked & (0x01 << i)) {
            this->outputs[i].appendBit(1);
        } else {
            this->outputs[i].appendBit(0);
        }
    }

    // Add new pattern to pattern map
    this->patternMap.insert(make_pair(iPatternMasked, oPatternMasked));
}


// Add a pattern as a pair
void truthTable::addPattern(pair<uint32_t, uint32_t> pattern) {
    this->addPattern(pattern.first, pattern.second);
}


// Get a pattern from the
pair<uint32_t, uint32_t> truthTable::getPattern(uint32_t index) {
    uint32_t inputBitmap = 0;
    uint32_t outputBitmap = 0;

    // Get input bitmap
    for(unsigned i = 0; i < this->getInputCount(); i++) {
        if(this->inputs[i].getBit(index)) {
            inputBitmap |= 0x01 << i;
        }
    }

    // Get output bitmap
    for(unsigned i = 0; i < this->getOutputCount(); i++) {
        if(this->outputs[i].getBit(index)) {
            outputBitmap |= 0x01 << i;
        }
    }

    // Make the pair and return
    return make_pair(inputBitmap, outputBitmap);
}


// Gets the input bitmap associated with
uint64_t truthTable::getInputBitmap(uint32_t inputIndex, uint32_t bitmapIndex) {
    return this->inputs[inputIndex].getBitmap(bitmapIndex);
}


// Gets the output bitmap associated with
uint64_t truthTable::getOutputBitmap(uint32_t outputIndex, uint32_t bitmapIndex) {
    return this->outputs[outputIndex].getBitmap(bitmapIndex);
}


// Gets bitmap mask for end bitmaps
uint64_t truthTable::getBitmapMask(uint32_t bitmapIndex) {
    return this->inputs[0].bitmapMask(bitmapIndex);
}


// Writes the truth table to a file with the specified radix
void truthTable::writeToFile(string path, uint32_t radix) {

    // Open file for writing
    ofstream fp(path);

    // Check that the file is indeed open
    if(!fp.is_open()) {
        throw(runtime_error("Could not open file '" + path + "' for writing."));
    }

    // Write radix and input/output count
    fp << "radix " << 2 << ";\n";
    fp << "iCount " << this->getInputCount() << ";\n";
    fp << "oCount " << this->getOutputCount() << ";\n";

    // Write patterns to file
    for(unsigned i = 0; i < this->getPatternCount(); i++) {
        pair<uint32_t, uint32_t> pattern = this->getPattern(i);
        fp << "pattern ";
        for(int j = this->getInputCount() - 1; j > -1; j--) {
            if(pattern.first & (0x01 << j)) fp << "1"; else fp << "0";
        }
        fp << ":";
        for(int j = this->getOutputCount() - 1; j > -1; j--) {
            if(pattern.second & (0x01 << j)) fp << "1"; else fp << "0";
        }
        fp << ";\n";
    }

    // Close the file
    fp.close();
}


// Write to file with default parameters
void truthTable::writeToFile(string path) {
    this->writeToFile(path, 2);
}
