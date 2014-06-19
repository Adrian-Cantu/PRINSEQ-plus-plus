//
//  Fasta.h
//  PrinSeq
//
//  Created by Jeffrey Sadural on 5/28/14.
//  Copyright (c) 2014 Jeffrey Sadural. All rights reserved.
//

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
#include "FormatCheck.h"
#include "SequenceData.h"

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
    
    ///Output Functions
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
    void TrimTailLeft();
    void TrimTailRight();
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
    
    string inputFileName;
    string badFileName;
    string goodFileName;
    ifstream fastaFile;
    ofstream GoodFileStream;
    ofstream BadFileStream;
    
    FormatCheck CheckFileFormat;
    SequenceData FastaSeq;
    
};



#endif /* defined(__PrinSeq__Fasta__) */
