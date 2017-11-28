#ifndef IOSTREAM
#define IOSTREAM
#include <iostream>
#endif

#ifndef REGEX
#define REGEX
#include <regex>
#endif

#include "reads.h"

using namespace std;
        single_read::single_read(istream &is): file1(is)  { 
            fastq_to_fasta.assign("^@");
        }

        void single_read::set_outputs(ostream& bad_out_file, ostream& single_out_file, ostream& good_out_file) {
            bad_out=bad_out_file.rdbuf();
            single_out=single_out_file.rdbuf();
            good_out=good_out_file.rdbuf();
        }

        int single_read::read_read(void) {
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

        int single_read::ns_max_n(int ns_max_n) {
            int hit_num=0;
            for( std::string::size_type i = seq_qual.size(); i > 0; --i) {
                if ((seq_seq[i-1] == 'n') || (seq_seq[i-1] == 'N')) {
                    hit_num++;
                }    
            }
            if ( hit_num > ns_max_n ) {
              read_status=2;
                return 1;
            } else {
                return 0;
            }
        }

        void single_read::print(int out_form) {
            if (read_status==2) {
                cout.rdbuf(bad_out);
            } else if (read_status==1) {
                cout.rdbuf(single_out);
            } else if (read_status==0) {
                cout.rdbuf(good_out);
            }
            if (out_form==0) {
                cout << seq_name << endl << seq_seq << endl << seq_sep << endl << seq_qual << endl;
            } else if (out_form==1) {
                cout << regex_replace(seq_name,fastq_to_fasta, ">") << endl << seq_seq << endl;
            }
            cout.rdbuf(back_stdout);
        }

        int single_read::min_qual_score(int min_qual) {
            string temp_seq_qual=seq_qual;
            int score;
            int i;
            for(i = seq_qual.size()-1; i >= 0; --i) {
                score=int(seq_qual[i])-33;
                if (score < min_qual) { return 1 ;}
            }
            return 0;
        }    

        int single_read::get_read_status(void) {
            return read_status;
        }

        void single_read::set_read_status(int status) {
            if (status > read_status) {read_status=status;}
        }

int single_read::min_qual_mean(int min_qual) {
    int score;
    float average=0;
    for(std::string::size_type i = seq_qual.size()-1; i > 0; --i) {
        score=int(seq_qual[i])-33;
        average= average + score;
        }
    average=average/seq_qual.size();
    if (average < min_qual) { return 1 ;}
    return 0;
}    

int single_read::noiupac() {
    regex pattern("^[ACGTN]+$", regex::icase);
    if (!regex_search(seq_seq,pattern)) {
        return 1;
    } else {
        return 0;
    }
}    


int single_read::min_len(unsigned int len) {
    if (seq_seq.size() < len) {
        return 1;
    } else {
        return 0;
    }
}

int single_read::max_len(unsigned int len) {
    if (seq_seq.size() > len) {
        return 1;
    } else {
        return 0;
    }
}





int single_read::max_gc(float max_gc) {
    int hit_num=0;
    for( std::string::size_type i = seq_qual.size(); i > 0; --i) {
        if ((seq_seq[i-1] == 'G') || (seq_seq[i-1] == 'C')
         || (seq_seq[i-1] == 'g') || (seq_seq[i-1] == 'c')) {
            hit_num++;
        }    
    }
    cout << "max_gc: " << max_gc << " , percent : " <<
    hit_num <<" / " << seq_seq.size() << " = " <<
    100*(float)hit_num/seq_seq.size()<< endl;
    if (max_gc < 100*(float)hit_num/seq_seq.size()) { return 1 ;}
    return 0;
}


int single_read::min_gc(float min_gc) {
    int hit_num=0;
    for( std::string::size_type i = seq_qual.size(); i > 0; --i) {
        if ((seq_seq[i-1] == 'G') || (seq_seq[i-1] == 'C')
         || (seq_seq[i-1] == 'g') || (seq_seq[i-1] == 'c')) {
            hit_num++;
        }    
    }
    if (min_gc > 100*(float)hit_num/seq_seq.size()) { return 1 ;}
    return 0;
}


//////////////////////////////////////////////////////////////////////////////

        pair_read::pair_read(istream &is1, istream &is2): file1(is1),file2(is2)  {

        read1= new single_read(file1);
        read2= new single_read(file2);
    }

    int pair_read::read_read(void) {
        return read1->read_read() * read2->read_read();
    }


    void pair_read::print(void) {
        read1->print(out_form);
        read2->print(out_form);
    }


    void pair_read::set_outputs(ostream& bad_out_file1, ostream& single_out_file1, ostream& good_out_file1,
                    ostream& bad_out_file2, ostream& single_out_file2, ostream& good_out_file2) {
        read1->set_outputs(bad_out_file1, single_out_file1, good_out_file1);
        read2->set_outputs(bad_out_file2, single_out_file2, good_out_file2);
    }

    void pair_read::ns_max_n(int ns_max_n) {
        regex pattern("n", regex::icase);
        int match1= read1->ns_max_n(ns_max_n);
        int match2= read2->ns_max_n(ns_max_n);
        pair_read::set_read_status(match1,match2);
    }

    void pair_read::min_qual_score(int min_qual) {
        int match1= read1->min_qual_score(min_qual);
        int match2= read2->min_qual_score(min_qual);
        pair_read::set_read_status(match1,match2);
    }

    void pair_read::min_qual_mean(int min_qual) {
        int match1= read1->min_qual_mean(min_qual);
        int match2= read2->min_qual_mean(min_qual);
        
        pair_read::set_read_status(match1,match2);
  }

void pair_read::noiupac(void) {
    int match1= read1->noiupac();
    int match2= read2->noiupac();
    pair_read::set_read_status(match1,match2);
}    


void pair_read::min_len(unsigned int len) {
    int match1= read1->min_len(len);
    int match2= read2->min_len(len);
    pair_read::set_read_status(match1,match2);
}    


void pair_read::max_len(unsigned int len) {
    int match1= read1->max_len(len);
    int match2= read2->max_len(len);
    pair_read::set_read_status(match1,match2);
}    


void pair_read::max_gc(float max_gc) {
    int match1= read1->max_gc(max_gc);
    int match2= read2->max_gc(max_gc);
    pair_read::set_read_status(match1,match2);
}

void pair_read::min_gc(float min_gc) {
    int match1= read1->min_gc(min_gc);
    int match2= read2->min_gc(min_gc);
    pair_read::set_read_status(match1,match2);
}
    void pair_read::set_out_format(int format) {
        out_form=format;
    }    

    void pair_read::set_read_status(int match1, int match2) {
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



string random_string( size_t length )
{
    auto randchar = []() -> char
    {
        const char charset[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
        const size_t max_index = (sizeof(charset) - 1);
        return charset[ rand() % max_index ];
    };
    std::string str(length,0);
    srand (time(NULL) ^ (time_t)&str );
    std::generate_n( str.begin(), length, randchar );
    return str;
}
