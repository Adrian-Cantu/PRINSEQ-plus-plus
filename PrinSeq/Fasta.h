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
#include <boost/regex.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>
#include "FormatCheck.h"

using std::string;
using std::fstream;
using std::ifstream;
using std::vector;
namespace po = boost::program_options;
using namespace std;

class Fasta{
public:
    Fasta();
    Fasta(int argc, char *fileLoc[]);
    void ParseOptions(int count, char *fileLoc[]);
    
    string RandFN();
    int DefaultOuputType(string name);
    void WriteGood(string filename, bool amino);
    void WriteBad(string filename, bool amino);
   
    //int CountFASeq(char fileLoc[], bool amino);
    //int CountFABase();
    //float meanLength();
    //void SetAmino(bool amino);
    //void WriteFasta();
    //bool ValExp(string s);
    
private:
    bool aa;  //amino acid flag
    int* base;
    int seqCount;
    int numbase;
    int outFormat;
    int argc;
    //char *fileLoc[];
    string fileName;
    string bOutFN;
    string gOutFN;
    
    FormatCheck Cff;
};

#endif /* defined(__PrinSeq__Fasta__) */
