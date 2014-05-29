//
//  FormatCheck.h
//  PrinSeq
//
//  Created by Jeffrey Sadural on 5/29/14.
//  Copyright (c) 2014 Jeffrey Sadural. All rights reserved.
//

#ifndef __PrinSeq__FormatCheck__
#define __PrinSeq__FormatCheck__

#include <iostream>
#include<fstream>
#include <boost/regex.hpp>
//#include<vector>
using std::string;
//using std::vector;
using std::fstream;
using std::ifstream;
class FormatCheck{
	
public:
	FormatCheck();
	string CheckFormat(string filename, bool amino);
	//int CountFASeq();
	//int CountFABase();
	//int CountQualSeq();
	//int CountQualBase();
	//int CountFQSeq();
	//int CountFQBase();
	//float meanLength();
	void SetAmino(bool amino);
	//bool FillList(char fileLoc[]);
	//bool FillFAList(char fileLoc[]);
	//bool FillQual(char fileLoc[]);
	//void WriteFasta();
	//void WriteQual();
	//void WriteFastq();
	
private:
	int fasta;
	int fastq;
	int qual;
	//vector<string>::iterator i;
	//vector<string> qualList;
	//vector<string> seqList;
	char name[256];
	int opt;
	string format;
	string s;
	string lineStore;
	bool aa;  //amino acid flag
	int found; //controls invalid character finder
	ifstream indata;
};
#endif /* defined(__PrinSeq__FormatCheck__) */
