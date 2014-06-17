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
    
    /// Option functions
    void DefineOptions(int optionCount, char *OptionsArray[]);
    void ProcessOptions();
    
    /// Not sure
    string RandFN();
    void ProcessData();
    bool ValExp(string s);
    
    ///Output Functions
    int DefaultOuputType(string name);
    void SetOutputFormat(int format);
    void WriteToGood(string filename);
    void WriteToBad(string filename);
    
    /// Summary Statistics
    void AddSeqCount();
    void AddBaseCount(int size);
    int GetBaseCount();
    int GetSeqCount();
    void ProcessFile();
    void Stats_All();
    
    /// Filters
    void ApplyFilters();
    void MinLengthFilter(string name);
    void MaxLengthFilter(string name);
    
    
    
   
    //int CountFASeq(char fileLoc[], bool amino);
    //int CountFABase();
    //float meanLength();
    //void SetAmino(bool amino);
    //void WriteFasta();
    //bool ValExp(string s);
    
private:
    
    po::options_description desc();
    
    bool amino;  //amino acid flag
    //int* base;
    int seqCount;
    int baseCount;
    int outFormat;
    //int argc;
    int seqNum;
    int trimLeft;
    int trimRight;
    int trimQualLeft;
    int trimQualRight;
    int trimTailLeft;
    int trimTailRight;
    int trimNSLeft;
    int trimNSRight;
    int trimToLen;
    int minLength;
    int maxLength;
    
    bool good;
    po::variables_map vm;
    
    //char *fileLoc[];
    string fileName;
    string badFileName;
    string goodFileName;
    ifstream fastaFile;
    ofstream GoodFileStream;
    ofstream BadFileStream;
    
    FormatCheck CheckFileFormat;
    SequenceData FastaSeq;
    
};



#endif /* defined(__PrinSeq__Fasta__) */
