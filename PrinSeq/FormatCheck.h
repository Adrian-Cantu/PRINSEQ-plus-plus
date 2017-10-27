//
//  FormatCheck.h
//  PrinSeq
//
//  Created by Jeffrey Sadural on 5/29/14.
//  Copyright (c) 2014 Jeffrey Sadural. All rights reserved.
//

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

