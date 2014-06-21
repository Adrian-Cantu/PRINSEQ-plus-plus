//
//  Options.h
//  PrinSeq
//
//  Created by Jeffrey Sadural on 6/21/14.
//  Copyright (c) 2014 Jeffrey Sadural. All rights reserved.
//

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
    
private:
    po::variables_map vm;
    
};

#endif /* defined(__PrinSeq__Options__) */
