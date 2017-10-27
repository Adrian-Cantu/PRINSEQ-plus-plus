//
//  FormatCheck.cpp
//  PrinSeq
//
//  Created by Jeffrey Sadural on 5/29/14.
//  Copyright (c) 2014 Jeffrey Sadural. All rights reserved.
//


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
