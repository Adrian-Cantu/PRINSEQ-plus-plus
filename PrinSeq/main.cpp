//===============================================================================
//   Authors: Robert SCHMIEDER and Jeff SADURAL
//	 Computational Science Research Center @ SDSU, CA
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

#include <boost/program_options.hpp>
#include "Fasta.h"
#include "Qual.h"
#include "SequenceData.h"
namespace po = boost::program_options;

#include <iostream>
#include <iterator>
#include <regex>
#include <string>
using namespace std;

int main(int argc, char* argv[]){
    
    srand((unsigned)time(0));
	int optind = 1;
    //Fasta FA(argc, argv);
	if (strcmp(argv[optind],"-fasta")){
        if (argv[optind + 1] == NULL) {
            cout << "Option fasta requires an argument." << endl;
            return 0;
        }
        
        if (argv[optind + 2] == NULL) {
            cout << "Nothing to do with input data." << endl;
            return 0;
        }
		Fasta FA(argc, argv);
        FA.ProcessFile();
	}
    
    else
        cout << "Error: you did not specify an input file containing the query sequences." << endl;
    
    return 0;
}
