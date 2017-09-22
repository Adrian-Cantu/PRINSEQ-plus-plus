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
  int index;
  int c;

  opterr = 0;


  while ((c = getopt (argc, argv, "abf:r:")) != -1)
    switch (c)
      {
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
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort ();
      }


  printf ("aflag = %d, bflag = %d, forward_read_file = %s\n",
          aflag, bflag, forward_read_file);

  for (index = optind; index < argc; index++)
    printf ("Non-option argument %s\n", argv[index]);

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
  
  regex pattern("n", regex::icase);
  ofstream bad_out_file;
  ofstream single_out_file;
  ofstream good_out_file;
  bad_out_file.open("test_bad_out.fastq");
  single_out_file.open("test_single_out.fastq");
  good_out_file.open("test_good_out.fastq");

  single_read read_f(inFile_f);
  single_read read_r(inFile_r);
  read_f.set_outputs(bad_out_file,single_out_file,good_out_file);
  read_r.set_outputs(bad_out_file,single_out_file,good_out_file);

  while(read_f.read_read()) {
    read_f.seq_match(pattern);
    read_f.print();
  }  
  while(read_r.read_read()) {
    read_r.seq_match(pattern);
    read_r.print();
  }
inFile_f.close();  
inFile_r.close();  
  return 0;
}
