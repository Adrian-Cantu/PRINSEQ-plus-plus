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

        int get_read_status(void) {
            return read_status;
        }

        void set_read_status(int status) {
            read_status=status;
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

class pair_read {
    public:
        pair_read(istream &is1, istream &is2): file1(is1),file2(is2)  {

        read1= new single_read(file1);
        read2= new single_read(file2);
    }

    int read_read(void) {
        return read1->read_read() * read2->read_read();
    }


    void print(void) {
        read1->print();
        read2->print();
    }


    int set_outputs(ostream& bad_out_file1, ostream& single_out_file1, ostream& good_out_file1,
                    ostream& bad_out_file2, ostream& single_out_file2, ostream& good_out_file2) {
        read1->set_outputs(bad_out_file1, single_out_file1, good_out_file1);
        read2->set_outputs(bad_out_file2, single_out_file2, good_out_file2);
    }

    int seq_match(regex pattern) {
        int match1= read1->seq_match(pattern);
        int match2= read2->seq_match(pattern);
        if ( !match1 && !match2 ) {
            read1->set_read_status(0);
            read2->set_read_status(0);
        } else if ( !match1 && match2 ) {
            read1->set_read_status(1);
            read2->set_read_status(2);
        } else if ( match1 && !match2 ) {
            read1->set_read_status(2);
            read2->set_read_status(1);
        } else { 
            read1->set_read_status(2);
            read2->set_read_status(2);
        }
    }



    protected:
        istream& file1;
        istream& file2;
        single_read* read1;
        single_read* read2;


};                    
