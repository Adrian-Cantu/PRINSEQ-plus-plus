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
        single_read(istream &is);
        int set_outputs(ostream& bad_out_file, ostream& single_out_file, ostream& good_out_file);
        int read_read(void);
        int seq_match(regex pattern,int ns_max_n);
        int max_n_p(int ns_max_p);
        void print(int out_form);
        int min_qual_score(int min_qual);
        int min_qual_mean(int min_qual);
        int get_read_status(void);
        void set_read_status(int status); 
        int noiupac(void);        


    protected:
        regex fastq_to_fasta;
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
        pair_read(istream &is1, istream &is2);
        int read_read(void);
        void print(void);
        int set_outputs(ostream& bad_out_file1, ostream& single_out_file1, ostream& good_out_file1,
                    ostream& bad_out_file2, ostream& single_out_file2, ostream& good_out_file2);
        int seq_match(regex pattern,int ns_max_n);
        int min_qual_score(int min_qual);
        int min_qual_mean(int min_qual);
        void set_out_format(int format);
        int max_n_p(int ns_max_p);
        void set_read_status(int match1, int match2);
        void noiupac(void);

    protected:
        istream& file1;
        istream& file2;
        single_read* read1;
        single_read* read2;
        int out_form=0;    //0 fastq , 1 fasta

};                    

string random_string( size_t length );
