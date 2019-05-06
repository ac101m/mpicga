#ifndef OPTPARSE_HPP_INCLUDED
#define OPTPARSE_HPP_INCLUDED


// Standard
#include <vector>
#include <string>
#include <sstream>


// Option type enum
typedef enum {
  ARG_TYPE_STRING,
  ARG_TYPE_INT,
  ARG_TYPE_FLOAT,
  ARG_TYPE_BOOL,
  ARG_TYPE_VOID
} ArgType;


// Conversion error exception
class ArgumentConversionException : std::exception {
  const char * what() const throw() {
    return "Conversion exception.";
  }
};


// Class to facilitate simple argument conversion
class Argument {
  private:

    // Base storage is a string
    std::string str;

  public:

    // Constructors
    Argument(std::string str);
    Argument(int i);
    Argument(bool b);
    Argument(double d);

    // Overloaded conversions
    std::string GetString();
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
    std::string longIdent;
    std::string desc;

    // Data retrieved from command line parameters
    std::vector<Argument> args;
    bool specified;

  public:

    // Constructors and initialiser
    Option(std::string longIdent, char shortIdent, ArgType type, std::string desc);
    Option(std::string longIdent, char shortIdent, ArgType type, std::string desc, std::vector<std::string> def);

    // Basic initialiser, for shared code
    void Init(std::string longIdent, char shortIdent, ArgType type, std::string desc);

    // Compare long and short identifiers
    bool OptMatches(char opt);
    bool OptMatches(std::string opt);

    // Parsing options
    void Specify(std::vector<std::string> args);

    // Get arguments in string form
    void AssertArgRequestValid(ArgType requestedType);
    bool Specified(void) {return this->specified;}

    // Vector conversions
    operator std::vector<std::string>();
    operator std::vector<int>();
    operator std::vector<bool>();
    operator std::vector<double>();

    // Single conversions
    operator std::string();
    operator int();
    operator bool();
    operator double();

    // Identifying string
    std::string IDStr(void);
    std::string ArgStr(void);
    std::string SyntaxStr(void);
    std::string HelpStr(void);
};


// Parser for command line arguments
class OptionParser {
  private:

    // Raw arguments, from command line
    std::vector<std::string> argv;

    // Options class and flags
    std::vector<Option> options;
    bool finalised;

    // Program description information
    std::string name;
    std::string desc;

  private:

    Option& FindOption(std::string longIdent);
    Option& FindOption(char shortIdent);

  public:

    // Constructors, and init routine
    void Init(int argc, char **argv);
    OptionParser(int argc, char **argv);
    OptionParser(int argc, char **argv, std::string desc);

    // Add an option to the parser
    void Add(Option opt);

    // Parsing routines
    void ParseLongOption(unsigned optIndex);
    void ParseShortOptionBlock(unsigned optIndex);
    void DoHelpOutput(void);
    void Parse(void);

    // Is a given option specified
    bool Specified(std::string longIdent);

    // Get a specified option
    Option Get(std::string longIdent);
};


#endif // OPTPARSE_HPP_INCLUDED
