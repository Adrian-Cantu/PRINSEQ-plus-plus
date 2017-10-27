//===============================================================================
//   Authors: Robert SCHMIEDER, Jeff SADURAL and Adrian CANTU
//     Computational Science Research Center @ SDSU, CA
//
//   File: prinseq-lite
//   Date: 2012-05-28
//   Version: 0.19.2 lite
//
//   Usage:
//      prinseq-lite [options]
//
//      Try 'prinseq-lite -h' for more information.
//
//   Purpose: PRINSEQ will help you to preprocess your genomic or metagenomic
//             sequence data in FASTA or FASTQ format. The lite version does not
//             require any non-core perl modules for processing.
//
//    Bugs: Please use http://sourceforge.net/tracker/?group_id=315449
//
//===============================================================================

#ifndef __PrinSeq__Options__
#define __PrinSeq__Options__

#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/regex.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <iomanip>

using std::string;
using std::fstream;
using std::ifstream;
using std::vector;
namespace po = boost::program_options;
using namespace std;

class Options{
public:
    Options();
    Options(int numberOfOptions, char *optionsArray[]);
    
    void DefineOptions(int numberOfOptions, char *OptionsArray[]);
    void ProcessOptions();
    bool IsOptionPresent(string option);
    string GetStringValue(string option);
    int GetIntValue(string option);
    bool GetBoolValue(string option);
    
private:
    po::variables_map vm;
    
};

#endif /* defined(__PrinSeq__Options__) */
