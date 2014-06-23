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


namespace po = boost::program_options;
using namespace std;

class Qual{
public:
    Qual();
    Qual(string label, string sequence);
    
    void ConvertQualNumsToAscii(string qualData);
    void ConvertToNumbers(string qualData);
    
private:
    
};
#endif /* defined(__PrinSeq__Qual__) */
