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
