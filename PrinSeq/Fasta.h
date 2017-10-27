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

#ifndef __PrinSeq__Fasta__
#define __PrinSeq__Fasta__

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
#include <exception>
#include "FormatCheck.h"
#include "SequenceData.h"
#include "Options.h"
#include "Qual.h"

using std::string;
using std::fstream;
using std::ifstream;
using std::vector;
namespace po = boost::program_options;
using namespace std;

class Fasta{
public:
    /// Constructors
    Fasta();
    Fasta(int optionCount, char *OptionsArray[]);
    void SetDefaultValues();
    
    /// Option functions
    void DefineOptions(int optionCount, char *OptionsArray[]);
    void ProcessOptions();
    
    /// Not sure
    string RandFN();
    void ProcessData();
    
    bool IsSeqID(string s);
    
    /// Output Functions
    int DefaultOuputType(string name);
    void SetOutputFormat(int format);
    void WriteToGood();
    void WriteToBad();
    
    /// Summary Statistics
    void IncrementSeqCount();
    void IncrementBaseCount(long size);
    void IncrementBadSeqCount();
    void IncrementGoodSeqCount();
    void IncrementBadBaseCount();
    void IncrementGoodBaseCount();
    
    long GetBaseCount();
    long GetSeqCount();
    long GetBadSeqCount();
    long GetGoodSeqCount();
    long GetBadBaseCount();
    long GetGoodBaseCount();
    
    void ProcessFile();
    
    /// Print Statistics
    void PrintStats();
    void PrintStandardStats();
    void PrintStatsInfo();
    void PrintStats_All();
    
    
    /// Filters
    void TrimSequence();
    void TrimQualLeft();
    void TrimQualRight();
    string CreateTail(char ATN, int tailLength);
    void TrimTailLeft(char base);
    void TrimTailRight(char base);
    void TrimNSLeft();
    void TrimNSRight();
    void ApplyFilters();
    void MinLengthFilter();
    void MaxLengthFilter();
    
private:
    
    po::options_description desc();
    
    bool amino;  //amino acid flag
    
    int outFormat;
    int seqNum;
    int trimLeftAmnt;
    int trimRightAmnt;
    int trimQualLeft;
    int trimQualRight;
    int trimTailLeft;
    int trimTailRight;
    int trimNSLeft;
    int trimNSRight;
    int trimToLen;
    int minLength;
    int maxLength;
    
    long seqCount;
    long baseCount;
    long goodSeqCount;
    long badSeqCount;
    long goodBaseCount;
    long badBaseCount;
    
    po::variables_map vm;
    Options optionMap;
    Qual qual;
    
    string inputFileName;
    string badFileName;
    string goodFileName;
    string tail;
    ifstream fastaFile;
    ofstream GoodFileStream;
    ofstream BadFileStream;
    
    FormatCheck CheckFileFormat;
    SequenceData FastaSeq;
    
};



#endif /* defined(__PrinSeq__Fasta__) */
