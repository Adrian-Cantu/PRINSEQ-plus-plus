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

#ifndef __PrinSeq__Qual__
#define __PrinSeq__Qual__

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
#include "Options.h"

namespace po = boost::program_options;
using namespace std;

class Qual{
public:
    Qual();
    Qual(string label, string sequence);
    
    string ConvertQualNumsToAscii(string qualData);
    int *ConvertASCIItoNums(string qualData);
    void ConvertToNumbers(string qualData);
    void trimFromLeft(int value);
    void trimFromRigh(int);
    
private:
    string qualSeq;
    int seqLength;
    vector<int> qualSeqArray;
    
    int trim_qual_left;
    int trim_qual_right;
    int trim_qual_window;
    int trim_qual_step;
    int qualityScore;
    string qualityType;
    string trim_qual_type;
    string trim_qual_rule;
    Options optionMap;
};
#endif /* defined(__PrinSeq__Qual__) */
