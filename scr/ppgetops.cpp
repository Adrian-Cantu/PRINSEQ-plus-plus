#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <fstream>


using namespace std;

int main (int argc, char **argv)
{
  int aflag = 0;
  int bflag = 0;
  char *cvalue = NULL;
  int index;
  int c;

  opterr = 0;


  while ((c = getopt (argc, argv, "abc:")) != -1)
    switch (c)
      {
      case 'a':
        aflag = 1;
        break;
      case 'b':
        bflag = 1;
        break;
      case 'c':
        cvalue = optarg;
        break;
      case '?':
        if (optopt == 'c')
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


  printf ("aflag = %d, bflag = %d, cvalue = %s\n",
          aflag, bflag, cvalue);

  for (index = optind; index < argc; index++)
    printf ("Non-option argument %s\n", argv[index]);

  ifstream inFile;
  inFile.open(cvalue);
  if (!inFile) {
    cerr << "Error: can not opem test file.txt" ;
    return 1;
    }
  string seq_name;
  string seq_seq;
  string seq_sep;
  string seq_qual;

  while (getline(inFile, seq_name, '\n')) {
    getline(inFile, seq_seq, '\n');
    getline(inFile, seq_sep, '\n');
    getline(inFile, seq_qual, '\n');
    cout << seq_name << endl << seq_seq << endl << seq_sep << endl << seq_qual << endl;
  } 


inFile.close();  
  return 0;
}

