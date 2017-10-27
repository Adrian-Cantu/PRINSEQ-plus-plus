#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#ifndef IOSTREAM
#define IOSTREAM
#include <iostream>
#endif

#include <fstream>

#ifndef REGEX
#define REGEX
#include <regex>
#endif

#include "reads.h"
using namespace std;

int main (int argc, char **argv)
{
    int aflag = 0;
    int bflag = 0;
    char *forward_read_file = NULL;
    char *reverse_read_file = NULL;
    string rand_string;
    int index;
    
    int c;

    opterr = 0;
    
    // Readin inout from the command line
    while ((c = getopt (argc, argv, "abf:r:")) != -1)
        switch (c) {
            case 'a':
                aflag = 1;
                break;
            case 'b':
                bflag = 1;
                break;
            case 'f':
                forward_read_file = optarg;
                break;
            case 'r':
                reverse_read_file = optarg;
                break;
            case '?':
                if (optopt == 'f')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
                return 1;
            default:
                abort ();
        }

    rand_string=random_string(6);    
    printf ("aflag = %d, bflag = %d, forward_read_file = %s ,reverse_read_file =%s\n ", aflag, bflag, forward_read_file, reverse_read_file);
//    printf ("random string = %s \n", rand_string);
    cout << "random string " << rand_string  << endl ;
    for (index = optind; index < argc; index++)
        printf ("Non-option argument %s\n", argv[index]);



/////////// open input files
    ifstream inFile_f;
    ifstream inFile_r;
    inFile_f.open(forward_read_file);
    if (!inFile_f) {
        cerr << "Error: can not open " << forward_read_file  << endl ;
        return 1;
    }
    inFile_r.open(reverse_read_file);
    if (!inFile_r) {
        cerr << "Error: can not open " << reverse_read_file  << endl ;
        return 1;
    }
      
////////// open and name oupput files
    ofstream bad_out_file_R1;
    ofstream single_out_file_R1;
    ofstream good_out_file_R1;
    ofstream bad_out_file_R2;
    ofstream single_out_file_R2;
    ofstream good_out_file_R2;

    bad_out_file_R1.open("Test_" + rand_string  + "_bad_out_R1.fastq");
    single_out_file_R1.open("Test_" + rand_string  + "_single_out_R1.fastq");
    good_out_file_R1.open("Test_" + rand_string  + "_good_out_R1.fastq");
    bad_out_file_R2.open("Test_" + rand_string  + "_bad_out_R2.fastq");
    single_out_file_R2.open("Test" + rand_string  + "__single_out_R2.fastq");
    good_out_file_R2.open("Test_"+ rand_string + "_good_out_R2.fastq");


    regex pattern("n", regex::icase);
    single_read read_f(inFile_f);
    single_read read_r(inFile_r);

    pair_read read_rf(inFile_f,inFile_r);

    read_rf.set_outputs(bad_out_file_R1,single_out_file_R1,good_out_file_R1,
        bad_out_file_R2,single_out_file_R2,good_out_file_R2);
/*
  while(read_f.read_read()) {
    read_f.seq_match(pattern);
    read_f.print();
  }  
  while(read_r.read_read()) {
    read_r.seq_match(pattern);
    read_r.print();
  }
*/  
    while(read_rf.read_read()) {
        read_rf.seq_match(pattern);
        read_rf.print();
    }

    inFile_f.close();  
    inFile_r.close();  
    return 0;
}
