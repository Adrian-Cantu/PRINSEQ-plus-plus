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

 //       int change_output(istream &good, istream &bad) {
            
        int set_outputs(ostream& bad_out_file, ostream& single_out_file, ostream& good_out_file) {
            bad_out=bad_out_file.rdbuf();
            single_out=single_out_file.rdbuf();
            good_out=good_out_file.rdbuf();
        }

        int read_read(void) {
        read_status=0;
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
                read_status=2;
                return 1;
            } else {
                return 0;
            }
        }

        void print(void) {
            if (read_status==2) {
                cout.rdbuf(bad_out);
            } else if (read_status==1) {
                cout.rdbuf(single_out);
            } else if (read_status==0) {
                cout.rdbuf(good_out);
            }
            cout << seq_name << endl << seq_seq << endl << seq_sep << endl << seq_qual << endl;
            cout.rdbuf(back_stdout);
        }    

    protected:
        int read_status=0; //0 good, 1 single ,2 bad
        istream& file1;
        string seq_name;
        string seq_seq;
        string seq_sep;
        string seq_qual;
        streambuf *back_stdout=cout.rdbuf();
        streambuf *bad_out=cout.rdbuf();
        streambuf *single_out=cout.rdbuf();
        streambuf *good_out=cout.rdbuf();
};         

