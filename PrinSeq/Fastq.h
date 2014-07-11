//
//  Fastq.h
//  PrinSeq
//
//  Created by Jeffrey Sadural on 7/9/14.
//  Copyright (c) 2014 Jeffrey Sadural. All rights reserved.
//

#ifndef __PrinSeq__Fastq__
#define __PrinSeq__Fastq__
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

class Fastq{
public:
    /// Constructors
    Fastq();
    Fastq(int optionCount, char *OptionsArray[]);
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
    ifstream fastqFile;
    ofstream GoodFileStream;
    ofstream BadFileStream;
    
    FormatCheck CheckFileFormat;
    SequenceData FastqSeq;
    
};

#endif /* defined(__PrinSeq__Fastq__) */
