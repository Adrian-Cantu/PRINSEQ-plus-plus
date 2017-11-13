#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

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
    char *forward_read_file = NULL;
    char *reverse_read_file = NULL;
    string rand_string;
    string out_ext="fastq";
    int index;
    int out_format=0; // 0 fastq, 1 fasta
    int min_qual_score=0;
    int ns_max_n=0;
    int min_qual_mean=0;
    int noiupac=0;
    int c;

    opterr = 0;

    struct option longopts[] = {
        { "fastq"           , required_argument , NULL     , 1 },
        { "fastq2"          , required_argument , NULL     , 2 },
        { "out_format"      , required_argument , NULL     , 3 },
        { "min_qual_score"  , required_argument , NULL     , 4 },
        { "ns_max_n"        , required_argument , NULL     , 5 },
        { "min_qual_mean"   , required_argument , NULL     , 6 },
        { "noiupac"         , no_argument       , &noiupac , 1 }, 
        {0,0,0,0}
    };    



    // Readin inout from the command line
    while ((c = getopt_long_only(argc, argv, "",longopts, NULL)) != -1)
        switch (c) {
            case 1:
                forward_read_file = optarg;
                break;
            case 2:
                reverse_read_file = optarg;
                break;
            case 3:
                out_format = atoi(optarg);
                break;
            case 4:
                min_qual_score = atoi(optarg);
                break;
            case 5:
                ns_max_n = atoi(optarg);
                break;
            case 6:
                min_qual_mean = atoi(optarg);
                break;
            case 0:
                // getopt set a variable
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
    printf ("forward_read_file = %s ,reverse_read_file =%s\n ", forward_read_file, reverse_read_file);
    cout << "random string " << rand_string << " out format " << out_format  << endl ;
    cout << "ns_max_n " << ns_max_n << endl ;
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

    if (out_format == 1 ) { out_ext = "fasta";}
    bad_out_file_R1.open("Test_" + rand_string  + "_bad_out_R1." + out_ext );
    single_out_file_R1.open("Test_" + rand_string  + "_single_out_R1." + out_ext  );
    good_out_file_R1.open("Test_" + rand_string  + "_good_out_R1." + out_ext);
    bad_out_file_R2.open("Test_" + rand_string  + "_bad_out_R2." + out_ext );
    single_out_file_R2.open("Test_" + rand_string  + "_single_out_R2." + out_ext);
    good_out_file_R2.open("Test_"+ rand_string + "_good_out_R2."  + out_ext);


    single_read read_f(inFile_f);
    single_read read_r(inFile_r);

    pair_read read_rf(inFile_f,inFile_r);

    read_rf.set_outputs(bad_out_file_R1,single_out_file_R1,good_out_file_R1,
        bad_out_file_R2,single_out_file_R2,good_out_file_R2);
    
    read_rf.set_out_format(out_format);

    // main loop
    while(read_rf.read_read()) {
        if (ns_max_n > -1 ) {read_rf.ns_max_n(ns_max_n);}
        if (min_qual_mean)  {read_rf.min_qual_mean(min_qual_mean);}
        if (min_qual_score) { read_rf.min_qual_score(min_qual_score);}
        if (noiupac) {read_rf.noiupac();}
        read_rf.print();
    }

    inFile_f.close();  
    inFile_r.close();  
    return 0;
}



