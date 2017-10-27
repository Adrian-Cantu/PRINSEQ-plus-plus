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


#include <iostream>
#include <fstream>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>


#define VALID_BASES_AA "[ABCDEFGHIKLMNOPQRSTUVWYZXabcdefghiklmmopqrstuvwyzx*-]+"

#define VALID_BASES_NON_AA "[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+"


using namespace boost;
using std::string;
using std::fstream;
using std::ifstream;

class FormatCheck{
	
public:
	FormatCheck();
	string CheckFormat(string filename, bool amino);
	void SetAmino(bool amino);

private:
	int fasta;
	int fastq;
	int qual;

	char name[256];
	int opt;
	string format;
	string s;
	string lineStore;
	bool aa;  //amino acid flag
	int found; //controls invalid character finder
	ifstream indata;
};

