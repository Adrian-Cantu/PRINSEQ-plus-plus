#ifndef IOSTREAM
#define IOSTREAM
#include <iostream>
#endif

#ifndef REGEX
#define REGEX
#include <regex>
#endif

#ifndef MATH
#define MATH
#include <math.h>
#endif

#include "reads.h"
#include <unordered_map>
#include <algorithm> 


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

        void single_read::ns_max_n(int ns_max_n) {
            int hit_num=0;
            for( std::string::size_type i = seq_qual.size(); i > 0; --i) {
                if ((seq_seq[i-1] == 'n') || (seq_seq[i-1] == 'N')) {
                    hit_num++;
                }    
            }
            if ( hit_num > ns_max_n ) {
            single_read::set_read_status(2);
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
                //cout << regex_replace(seq_name,fastq_to_fasta, '>') << endl << seq_seq << endl;
                string seq_name_copy=seq_name;
                seq_name_copy[0]='>';
                cout << seq_name_copy << endl << seq_seq << endl << seq_sep << endl << seq_qual << endl;
            }
            cout.rdbuf(back_stdout);
        }

        void single_read::min_qual_score(int min_qual) {
            string temp_seq_qual=seq_qual;
            int score;
            int i;
            for(i = seq_qual.size()-1; i >= 0; --i) {
                score=int(seq_qual[i])-33;
                if (score < min_qual) { single_read::set_read_status(2);}
            }
        }    

        int single_read::get_read_status(void) {
            return read_status;
        }

        void single_read::set_read_status(int status) {
            if (status > read_status) {read_status=status;}
        }

void single_read::min_qual_mean(int min_qual) {
    int score;
    float average=0;
    for(std::string::size_type i = seq_qual.size()-1; i > 0; --i) {
        score=int(seq_qual[i])-33;
        average= average + score;
        }
    average=average/seq_qual.size();
    if (average < min_qual) { single_read::set_read_status(2);}
}    

void single_read::noiupac() {
    regex pattern("^[ACGTN]+$", regex::icase);
    if (!regex_search(seq_seq,pattern)) {
        single_read::set_read_status(2);
    }
}    


void single_read::min_len(unsigned int len) {
    if (seq_seq.size() < len) {
        single_read::set_read_status(2);
    }    
}

void single_read::max_len(unsigned int len) {
    if (seq_seq.size() > len) {
        single_read::set_read_status(2);
    } 
}

void single_read::max_gc(float max_gc) {
    int hit_num=0;
    for( std::string::size_type i = seq_qual.size(); i > 0; --i) {
        if ((seq_seq[i-1] == 'G') || (seq_seq[i-1] == 'C')
         || (seq_seq[i-1] == 'g') || (seq_seq[i-1] == 'c')) {
            hit_num++;
        }    
    }
//    cout << "max_gc: " << max_gc << " , percent : " <<
//    hit_num <<" / " << seq_seq.size() << " = " <<
//    100*(float)hit_num/seq_seq.size()<< endl;
    if (max_gc < 100*(float)hit_num/seq_seq.size()) { single_read::set_read_status(2);}
}

void single_read::min_gc(float min_gc) {
    int hit_num=0;
    for( std::string::size_type i = seq_qual.size(); i > 0; --i) {
        if ((seq_seq[i-1] == 'G') || (seq_seq[i-1] == 'C')
         || (seq_seq[i-1] == 'g') || (seq_seq[i-1] == 'c')) {
            hit_num++;
        }    
    }
    if (min_gc > 100*(float)hit_num/seq_seq.size()) { single_read::set_read_status(2);}
}

void single_read::entropy(float threshold) {
    unsigned int j=0;
    std::string window;
    vector<float> vals;
    while ( 1 ) {
        try {
            window = seq_seq.substr(j*32,64);
        } catch (const std::exception& e) {
            break;
        }
        if (window.size() < 15 ) {
            if (vals.size() == 0) {
                vals.push_back(0.0);
                break;
            } else {
                break;
            }
        }
        unordered_map<string, int> hashtable;
        for ( std::string::size_type i = 0; i < window.size()-2 ; i++  ) {
            string sub = window.substr(i,3);
            if (hashtable.count(sub)) {
                hashtable[sub]++;
            } else {
                hashtable.emplace(sub,1);
            }
        }
        double entropy=0;
        double l=window.size()-2;
        double k=min(64.0,l);
        for (auto &itr : hashtable) {
            entropy -= ((float)itr.second/l)*(log((float)itr.second/l)/log(k));
        }
        vals.push_back(entropy);
        j++;    
        
    }
    double mean = 1.0 * std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    if (mean < threshold ) {
        single_read::set_read_status(2);
    }
}

void single_read::dust(float threshold) {
    unsigned int j=0;
    std::string window;
    vector<float> vals;
//    std::cout << seq_seq << std::endl;
    while ( 1 ) {
        try {
            window = seq_seq.substr(j*32,64);
        } catch (const std::exception& e) {
            break;
        }
        if (window.size() < 15 ) {
            if (vals.size() == 0) {
                vals.push_back(62.0);
                break;
            } else {
                break;
            }
        }
        unordered_map<string, int> hashtable;
        
        for ( std::string::size_type i = 0; i < window.size()-2 ; i++  ) {
            string sub = window.substr(i,3);
            if (hashtable.count(sub)) {
                hashtable[sub]++;
            } else {
                hashtable.emplace(sub,1);
            }
        }
        double dust=0;
        double l=window.size()-2;
        for (auto &itr : hashtable) {
            dust += (((float)itr.second)*((float)itr.second-1))/(l-1);
        }
        vals.push_back(dust);
 //       cout << "dust " << j+1 << " is : " << (dust*0.5)/(31) << endl;
        j++;    
        
    }
    double mean = 1.0 * std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    mean = (mean*0.5)/(31);
//    cout << "total entropy is " << mean << endl ;
    if (mean > threshold ) {
        single_read::set_read_status(2);
    }    
}

 // type min* mean max sum // rule lt* gt eq 
void single_read::trim_qual_right(string type, string rule, int step, int window_size, float threshold ) {
    string window;
    string copy_seq=seq_seq;
    string copy_qual=seq_qual;
//    std::cout << seq_seq << std::endl;
    while ( 1 ) {
        try {
            window = copy_qual.substr(copy_qual.size()-window_size,window_size);
        } catch (const std::exception& e) {
            break;
        }
        int score;
        vector<float> vals;
        for(int i = window.size()-1; i >= 0; --i) {
                score=int(window[i])-33;
                vals.push_back(score);
        }
        float compare;
        if (type == "min") {
            compare = *min_element(vals.begin(), vals.end());
        } else if (type == "max") {
            compare = *max_element(vals.begin(), vals.end());
        } else if (type == "mean" ) {
            compare = accumulate( vals.begin(), vals.end(), 0.0)/vals.size();
        } else {
            compare = accumulate( vals.begin(), vals.end(), 0.0);
        }
//        std::cout << window << " has and average score of " << compare << "and threshold of " << threshold << std::endl;
        if ((rule == "lt") && (compare < threshold)) {
            size_t b_win = max(0,(int)copy_qual.size()-step);
            copy_qual.erase(b_win,copy_qual.size());
            copy_seq.erase(b_win,copy_qual.size());
//            std::cout << copy_seq << std::endl;
        } else {
            break;
        }
    }
    if (copy_qual.size() == 0) {
        single_read::set_read_status(2);
    } else {
        seq_qual=copy_qual;
        seq_seq=copy_seq;
    }
//    std::cout << seq_seq << std::endl << std::endl;
}


 // type min* mean max sum // rule lt* gt eq 
void single_read::trim_qual_left(string type, string rule, int step, int window_size, float threshold ) {
    string window;
    string copy_seq=seq_seq;
    string copy_qual=seq_qual;
//    std::cout << seq_seq << std::endl;
    while ( 1 ) {
        try {
            window = copy_qual.substr(0,window_size);
        } catch (const std::exception& e) {
            break;
        }
        int score;
        vector<float> vals;
        for(int i = window.size()-1; i >= 0; --i) {
                score=int(window[i])-33;
                vals.push_back(score);
        }
        float compare;
        if (type == "min") {
            compare = *min_element(vals.begin(), vals.end());
        } else if (type == "max") {
            compare = *max_element(vals.begin(), vals.end());
        } else if (type == "mean" ) {
            compare = accumulate( vals.begin(), vals.end(), 0.0)/vals.size();
        } else {
            compare = accumulate( vals.begin(), vals.end(), 0.0);
        }
  //      std::cout << window << " has and average score of " << compare << " and threshold of " << threshold << std::endl;
        if ((rule == "lt") && (compare < threshold)) {
            size_t b_win = min(step,(int)copy_qual.size());
            copy_qual.erase(0,b_win);
            copy_seq.erase(0,b_win);
//            std::cout << copy_seq << std::endl;
        } else {
            break;
        }
    }
    if (copy_qual.size() == 0) {
        single_read::set_read_status(2);
    } else {
        seq_qual=copy_qual;
        seq_seq=copy_seq;
    }
//    std::cout << seq_seq << std::endl << std::endl;
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
        read1->ns_max_n(ns_max_n);
        read2->ns_max_n(ns_max_n);
        pair_read::auto_set_read_status();
    }

    void pair_read::min_qual_score(int min_qual) {
        read1->min_qual_score(min_qual);
        read2->min_qual_score(min_qual);
        pair_read::auto_set_read_status();
    }

    void pair_read::min_qual_mean(int min_qual) {
        read1->min_qual_mean(min_qual);
        read2->min_qual_mean(min_qual);
        pair_read::auto_set_read_status();
    }

    void pair_read::noiupac(void) {
        read1->noiupac();
        read2->noiupac();
        pair_read::auto_set_read_status();
    }    


void pair_read::min_len(unsigned int len) {
    read1->min_len(len);
    read2->min_len(len);
    pair_read::auto_set_read_status();
}    


void pair_read::max_len(unsigned int len) {
    read1->max_len(len);
    read2->max_len(len);
    pair_read::auto_set_read_status();
}    


void pair_read::max_gc(float max_gc) {
    read1->max_gc(max_gc);
    read2->max_gc(max_gc);
    pair_read::auto_set_read_status();
}

void pair_read::min_gc(float min_gc) {
    read1->min_gc(min_gc);
    read2->min_gc(min_gc);
    pair_read::auto_set_read_status();
}

    void pair_read::set_out_format(int format) {
        out_form=format;
    }

void pair_read::entropy(float threshold) {
    read1->entropy(threshold);
    read2->entropy(threshold);
    pair_read::auto_set_read_status();
}

void pair_read::dust(float threshold) {
    read1->dust(threshold);
    read2->dust(threshold);
    pair_read::auto_set_read_status();
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

    void pair_read::auto_set_read_status(void) {
        if ((read1->get_read_status()==0) && (read2->get_read_status()==2)) {
            read1->set_read_status(1);
        } else if ((read2->get_read_status()==0) && (read1->get_read_status()==2)) {
            read2->set_read_status(1);
        }
    }

void pair_read::trim_qual_right(string type, string rule, int step, int window_size, float threshold ) {
    read1->trim_qual_right(type, rule, step, window_size, threshold);
    read2->trim_qual_right(type, rule, step, window_size, threshold);
    pair_read::auto_set_read_status();
} 

void pair_read::trim_qual_left(string type, string rule, int step, int window_size, float threshold ) {
    read1->trim_qual_left(type, rule, step, window_size, threshold);
    read2->trim_qual_left(type, rule, step, window_size, threshold);
    pair_read::auto_set_read_status();
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
