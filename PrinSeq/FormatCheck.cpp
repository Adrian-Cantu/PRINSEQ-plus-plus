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



#include "FormatCheck.h"
#include<iostream>
#include<fstream>
using std::ifstream;

/*
 *  Default Constructor
 */
FormatCheck::FormatCheck(){
	fasta = 0;
	fastq = 0;
	qual = 0;
	
	format = "unknown";
	aa = true;  //amino acid flag
	found = 0; //controls invalid character finder
}

/*
 *  Checks if file is in a compatible format
 */
string FormatCheck::CheckFormat(string filename, bool amino){
	indata.open(filename);
	aa = amino;
    
	if(!indata) { // file couldn't be opened
		string error = "File could not be opened";
		return error;
	}
    
	for (int count = 3; count > 0; count--) {
		indata.getline(name, 256);
		s = name;

		while(s == "" && (indata.peek() != EOF)){  //skips blank lines
			indata.getline(name,256);
			s = name;
		}
		
		if(fasta == 0 && s[0] == '>') {
			fasta = 1;
			qual = 1;
		}
		
		else if(fasta == 1) {
			if (aa){
				found = (int)s.find_first_not_of(VALID_BASES_AA);
				if(found == -1){
					fasta = 2;
				}
				
                if(qual == 1){
                    //cout << "in qual if" << s[0] << endl;
                    if(isdigit(s[0])){
                        qual = 2;
                    }
                }
            }
			else if(!aa){
				found = (int)s.find_first_not_of(VALID_BASES_NON_AA);
				
				if(found == -1){
					fasta = 2;
				}
				
				if(qual == 1){
					//cout << "in qual if" << s[0] << endl;
					if(isdigit(s[0])){
						qual = 2;
					}
				}
			}
		}
		
		else if(fastq == 0 && s[0] == '@'){
			//id = 1
			fastq = 1;
		}
		
		else if(fastq == 1){
            
			if (aa){
				found = (int)s.find_first_not_of(VALID_BASES_AA);
				
				if(found == -1){
					fastq = 2;
				}
			}
			
			else if(!aa){
				found = (int)s.find_first_not_of(VALID_BASES_NON_AA);
				
				if(found == -1){
					fastq = 2;
				}
			}
		}
		
		else if(fastq == 2 && s[0] == '+'){
			fastq = 3;
		}
	}
	
	indata.close();
	
	if(qual == 2) {
		format = "qual";
	}
	
	else if(fasta == 2){
		format = "fasta";
	}
	
	else if(fastq == 3) {
		format = "fastq";
	}
	return format;
}
