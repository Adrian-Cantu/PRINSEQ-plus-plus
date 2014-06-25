//
//  Qual.h
//  PrinSeq
//
//  Created by Jeffrey Sadural on 6/9/14.
//  Copyright (c) 2014 Jeffrey Sadural. All rights reserved.
//

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
