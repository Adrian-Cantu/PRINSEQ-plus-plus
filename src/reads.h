#ifndef IOSTREAM
#define IOSTREAM
#include <iostream>
#endif

#ifndef REGEX
#define REGEX
#include <regex>
#endif


using namespace std;
class single_read {
    public:

        single_read(istream &is): file1(is)  {
        }

        int read_read(void) {
            if (getline(file1,seq_name, '\n')) {
                getline(file1, seq_seq, '\n');
                getline(file1, seq_sep, '\n');
                getline(file1, seq_qual, '\n');
                return 1;
            } else {
                return 0;
            }
        }

        int seq_match(regex pattern) {
            if (regex_search(seq_seq,pattern)) {
                return 1;
            } else {
                return 0;
            }
        }

        void print(void) {
            cout << seq_name << endl << seq_seq << endl << seq_sep << endl << seq_qual << endl;
        }    

    protected:
        istream& file1;
        string seq_name;
        string seq_seq;
        string seq_sep;
        string seq_qual;
};         

