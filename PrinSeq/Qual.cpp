//
//  Qual.cpp
//  PrinSeq
//
//  Created by Jeffrey Sadural on 6/9/14.
//  Copyright (c) 2014 Jeffrey Sadural. All rights reserved.
//

#include "Qual.h"


Qual::Qual(){
    trim_qual_left = 0;
    trim_qual_right = 0;
    trim_qual_window = 1;
    trim_qual_step = 1;
    trim_qual_type = "min";
    trim_qual_rule = "lt";
}

Qual::Qual(string label, string sequence){
    
}

string Qual::ConvertQualNumsToAscii(string qualData){
    string convertedData = string();
    int qualNumber;
    stringstream ss(qualData);

    while(ss >> qualNumber){
        convertedData += (char)qualNumber;
    }
    return convertedData;
}

int *ConvertASCIItoNums(string qualData){
    int *numarray;
    
    return numarray;
}

void Qual::ConvertToNumbers(string qualData){
    
}