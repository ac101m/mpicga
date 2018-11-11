#ifndef OPTPARSE_HPP_INCLUDED
#define OPTPARSE_HPP_INCLUDED


#include <vector>
#include <string>
#include <sstream>
using namespace std;


// Option type enum
typedef enum {
  ARG_TYPE_STRING,
  ARG_TYPE_INT,
  ARG_TYPE_FLOAT,
  ARG_TYPE_BOOL,
  ARG_TYPE_VOID
} ArgType;


// Conversion error exception
class ArgumentConversionException : exception {
  const char * what() const throw() {
    return "Conversion exception.";
  }
};


// Class to facilitate simple argument conversion
class Argument {
  private:

    // Base storage is a string
    string str;

  public:

    // Constructors
    Argument(string str);
    Argument(int i);
    Argument(bool b);
    Argument(double d);

    // Overloaded conversions
    string GetString();
    int GetInt();
    bool GetBool();
    double GetDouble();
};


// Class represents a possible command line parameter in the code
class Option {
  private:

    // Initial data
    ArgType type;
    char shortIdent;
    string longIdent;
    string desc;

    // Data retrieved from command line parameters
    vector<Argument> args;
    bool specified;

  public:

    // Constructors and initialiser
    Option(string longIdent, char shortIdent, ArgType type, string desc);
    Option(string longIdent, char shortIdent, ArgType type, string desc, vector<string> def);

    // Basic initialiser, for shared code
    void Init(string longIdent, char shortIdent, ArgType type, string desc);

    // Compare long and short identifiers
    bool OptMatches(char opt);
    bool OptMatches(string opt);

    // Parsing options
    void Specify(vector<string> args);

    // Get arguments in string form
    void AssertArgRequestValid(ArgType requestedType);
    bool Specified(void) {return this->specified;}

    // Vector conversions
    operator vector<string>();
    operator vector<int>();
    operator vector<bool>();
    operator vector<double>();

    // Single conversions
    operator string();
    operator int();
    operator bool();
    operator double();

    // Identifying string
    string IDStr(void);
    string ArgStr(void);
    string SyntaxStr(void);
    string HelpStr(void);
};


// Parser for command line arguments
class OptionParser {
  private:

    // Raw arguments, from command line
    vector<string> argv;

    // Options class and flags
    vector<Option> options;
    bool finalised;

    // Program description information
    string name;
    string desc;

  private:

    Option& FindOption(string longIdent);
    Option& FindOption(char shortIdent);

  public:

    // Constructors, and init routine
    void Init(int argc, char **argv);
    OptionParser(int argc, char **argv);
    OptionParser(int argc, char **argv, string desc);

    // Add an option to the parser
    void Add(Option opt);

    // Parsing routines
    void ParseLongOption(unsigned optIndex);
    void ParseShortOptionBlock(unsigned optIndex);
    void DoHelpOutput(void);
    void Parse(void);

    // Is a given option specified
    bool Specified(string longIdent);

    // Get a specified option
    Option Get(string longIdent);
};


#endif // OPTPARSE_HPP_INCLUDED
