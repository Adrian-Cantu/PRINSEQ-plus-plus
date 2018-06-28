/** \brief  
 * 
 */


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

#ifndef PTHREAD
#define PTHREAD
#include <pthread.h>
#endif

#include "reads.h"
#include <unordered_map>
#include <algorithm> 

#ifndef NUMERIC
#define NUMERIC
#include <numeric>
#endif


using namespace std;
single_read::single_read(istream &is, int mode): file1(is) , qual_mode(mode)  { 
    fastq_to_fasta.assign("^@");
    out_stream = new ostream(nullptr);
}
        
single_read::single_read(void) : file1(cin){ // starndar input
    fastq_to_fasta.assign("^@");
    out_stream = new ostream(nullptr);
    qual_mode=33;
}
        
/** \brief Set or change inmput stream.
 * 
 * Mainlly used after the default constructor which set the the input stream
 * to std::cin
 * 
 */
void  single_read::set_inputs(istream &is) {
    file1.rdbuf(is.rdbuf());
}

/** \brief Set output stream for good, single and bad read.
 * 
 * good read are ones that pass all filters, bad reads dont pass at least one filter.
 * single reads are those that pass filters but mate didn't. Only used for pair end reads.
 */ 
void single_read::set_outputs(ostream& bad_out_file, ostream& single_out_file, ostream& good_out_file) {
    bad_out=bad_out_file.rdbuf();
    single_out=single_out_file.rdbuf();
    good_out=good_out_file.rdbuf();
}

/** \brief get a read from the input stream and reset read status
 * 
 */ 
int single_read::read_read(pthread_mutex_t * read_mutex,int format) {
    read_status=0;
    string fasta_sequence="";
    string temp_fasta="";
    std::string token;
    int ff=0;
    seq_name="";
    seq_seq="";
    if (format==1) {
        pthread_mutex_lock(read_mutex);
        LOOP: if (getline(file1,temp_fasta, '>')) {
            if (temp_fasta.empty() && (!fasta_sequence.empty() || ff)) {  // thid deals with >> on some part of the identifier
                fasta_sequence.push_back('>');
                goto LOOP;
            }
            if (temp_fasta.empty() && fasta_sequence.empty()) { 
                ff++; 
                goto LOOP;
            } // del mainly with the first sequence id

            if (temp_fasta.back() != '\n') { // deals with > in the id other than the one at the start of the line
                fasta_sequence.append(temp_fasta);
                fasta_sequence.push_back('>');
                goto LOOP;
            } else {
                fasta_sequence.append(temp_fasta);
            }
            stringstream fasta_stream(fasta_sequence);
            getline(fasta_stream,seq_name);
            seq_name.insert(0,"@");
            while (getline(fasta_stream,temp_fasta)) {
                seq_seq.append(temp_fasta);
            }
            seq_sep='+';
            seq_qual= string(seq_seq.size(),'A');
            pthread_mutex_unlock(read_mutex);
            return 1;
        } else {
            pthread_mutex_unlock(read_mutex);
            return 0;
        }
        
    } else {
        pthread_mutex_lock(read_mutex);
        if (getline(file1,seq_name, '\n')) {
            getline(file1, seq_seq, '\n');
            getline(file1, seq_sep, '\n');
            getline(file1, seq_qual, '\n');
            pthread_mutex_unlock(read_mutex);
            return 1;
        } else {
            pthread_mutex_unlock(read_mutex);
            return 0;
        }
    }
}

/** \brief Filter out reads with more n's than \p ns_max_n
 * 
 */
int single_read::ns_max_n(int ns_max_n) {
    if (read_status==2) {return 0;}
    int hit_num=0;
    for( std::string::size_type i = seq_qual.size(); i > 0; --i) {
        if ((seq_seq[i-1] == 'n') || (seq_seq[i-1] == 'N')) {
            hit_num++;
        }    
    }
    if ( hit_num > ns_max_n ) {
        single_read::set_read_status(2);
        return 1;
    }
    return 0;
}

/** \brief Print read to the apropiated output stream
 * 
 * Output stream are set by single_read::set_outputs. \p out_form define the
 * output format, 0 for FASTQ and 1 for FASTA.
 * 
 */
void single_read::print(int out_form) {
    if (read_status==2) {
        out_stream->rdbuf(bad_out);
    } else if (read_status==1) {
        out_stream->rdbuf(single_out);
    } else if (read_status==0) {
        out_stream->rdbuf(good_out);
    }
    if (out_form==0) {
        *out_stream << seq_name << endl << seq_seq << endl << seq_sep << endl << seq_qual << endl;
    } else if (out_form==1) {
        string seq_name_copy=seq_name;
        seq_name_copy[0]='>';
        *out_stream << seq_name_copy << endl << seq_seq << endl;
    }
    //cout.rdbuf(back_stdout);
}

/** \brief filter out reads with at least one base quality below \p min_qual
 * 
 */
int single_read::min_qual_score(int min_qual) {
    string temp_seq_qual=seq_qual;
    int score;
    int i;
    if (read_status==2) {return 0;}
    for(i = seq_qual.size()-1; i >= 0; --i) {
        score=int(seq_qual[i])-qual_mode;
        if (score < min_qual) { 
            single_read::set_read_status(2);
            return 1;
        }
    }
    return 0;
}    


int single_read::get_read_status(void) {
    return read_status;
}

/** \brief Change the read status if \p status is worst
 * 
 * can ghange good to bad or single, and single to bad.
 * 
 */
void single_read::set_read_status(int status) {
    if (status > read_status) {read_status=status;}
}

/** \brief filter out reads with mean base quality below \p min_qual
 * 
 */
int single_read::min_qual_mean(int min_qual) {
    int score;
    float average=0;
    if (read_status==2) {return 0;}
    for(std::string::size_type i = seq_qual.size()-1; i > 0; --i) {
        score=int(seq_qual[i])-qual_mode;
        average= average + score;
        }
    average=average/seq_qual.size();
    if (average < min_qual) { 
        single_read::set_read_status(2);
        return 1;
    }
    return 0;
}    

/** \brief Filter out reads with iupac extended bases
 * 
 */
int single_read::noiupac() {
    if (read_status==2) {return 0;}
    regex pattern("^[ACGTN]+$", regex::icase);
    if (!regex_search(seq_seq,pattern)) {
        single_read::set_read_status(2);
        return 1;
    }
    return 0;
}    


/** \brief Filter out reads shorther than \p len
 * 
 */
int single_read::min_len(unsigned int len) {
    if (read_status==2) {return 0;}
    if (seq_seq.size() < len) {
        single_read::set_read_status(2);
        return 1;
    } else {
        return 0;
    }   
}

/** \brief Filter out reads longer than \p len
 * 
 */
int single_read::max_len(unsigned int len) {
    if (read_status==2) {return 0;}
    if (seq_seq.size() > len) {
        single_read::set_read_status(2);
        return 1;
    }
    return 0;
}

/** \brief Filter out reads with gc% higher than \p max_gc
 * 
 */
int single_read::max_gc(float max_gc) {
    int hit_num=0;
    if (read_status==2) {return 0;}
    for( std::string::size_type i = seq_qual.size(); i > 0; --i) {
        if ((seq_seq[i-1] == 'G') || (seq_seq[i-1] == 'C')
         || (seq_seq[i-1] == 'g') || (seq_seq[i-1] == 'c')) {
            hit_num++;
        }    
    }
    if (max_gc < 100*(float)hit_num/seq_seq.size()) { 
        single_read::set_read_status(2);
        return 1;
    }
    return 0;
}


/** \brief Filter out reads with gc% lower than \p min_gc
 * 
 */
int single_read::min_gc(float min_gc) {
    if (read_status==2) {return 0;}
    int hit_num=0;
    for( std::string::size_type i = seq_qual.size(); i > 0; --i) {
        if ((seq_seq[i-1] == 'G') || (seq_seq[i-1] == 'C')
         || (seq_seq[i-1] == 'g') || (seq_seq[i-1] == 'c')) {
            hit_num++;
        }    
    }
    if (min_gc > 100*(float)hit_num/seq_seq.size()) { 
        single_read::set_read_status(2);
        return 1;
    }
    return 0;
}

/** \brief Filter out reads with information lower than \p threshold
 * 
 * Uses the Shanon-Wiener entropy
 * 
 * \f[ CE=- \sum_{i=1}^{k}\left ( \frac{n_i}{l} \right )log_k\left ( \frac{n_i}{l} \right )  \f]
 */
int single_read::entropy(float threshold) {
    if (read_status==2) {return 0;}
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
        return 1;
    }
    return 0;
}

int single_read::dust(float threshold) {
    if (read_status==2) {return 0;}
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
        return 1;
    }
    return 0;
}

 // type min* mean max sum // rule lt* gt eq 
int single_read::trim_qual_right(string type, string rule, int step, int window_size, float threshold ) {
    if (read_status==2) {return 0;}
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
                score=int(window[i])-qual_mode;
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
        return 1;
    } else {
        seq_qual=copy_qual;
        seq_seq=copy_seq;
        return 0;
    }
//    std::cout << seq_seq << std::endl << std::endl;
}


 // type min* mean max sum // rule lt* gt eq 
int single_read::trim_qual_left(string type, string rule, int step, int window_size, float threshold ) {
    if (read_status==2) {return 0;}
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
                score=int(window[i])-qual_mode;
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
        return 1;
    } else {
        seq_qual=copy_qual;
        seq_seq=copy_seq;
        return 0;
    }
//    std::cout << seq_seq << std::endl << std::endl;
}

void single_read::rm_header(void) {
    seq_sep="+";
}

int single_read::trim_tail_left(int num) {
    if (read_status==2) {return 0;}
    int sum=0;
    int temp_size = seq_seq.size();
    for(int i = 0; i < temp_size; i++) {
        if ((seq_seq[i]=='A') || (seq_seq[i]=='T') || (seq_seq[i]=='a') || (seq_seq[i]=='t')) {
            sum++;
        } else {
            break;
        }
    }
    if (sum == temp_size ) {
        single_read::set_read_status(2);
        return 1;
    } else if (sum >= num) {
        seq_seq.erase(0,sum);
        seq_qual.erase(0,sum);
    }
    return 0;
}

int single_read::trim_tail_right(int num) {
    if (read_status==2) {return 0;}
    int sum=0;
    int temp_size = seq_seq.size();
    for(int i = temp_size -1; i >= 0; i--) {
        if ((seq_seq[i]=='A') || (seq_seq[i]=='T') || (seq_seq[i]=='a') || (seq_seq[i]=='t')) {
            sum++;
        } else {
            break;
        }
    }
    if (sum == temp_size ) {
        single_read::set_read_status(2);
        return 1;
    } else if (sum >= num) {
        seq_seq.erase(temp_size-sum,temp_size);
        seq_qual.erase(temp_size-sum,temp_size);
    }
    return 0;
}


//////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////

        pair_read::pair_read(istream &is1, istream &is2, int mode): file1(is1),file2(is2)  {

        read1= new single_read(file1,mode);
        read2= new single_read(file2,mode);
    }
    
    pair_read::pair_read(void):file1(cin),file2(cin) {
        read1= new single_read(file1,33);
        read2= new single_read(file2,33);
    }    
    
    void  pair_read::set_inputs(istream &read_f,istream &read_r) {
   //     read1->file1.rdbuf(read_f.rdbuf());
    //    read2->file1.rdbuf(read_r.rdbuf());
        read1->set_inputs(read_f);
        read2->set_inputs(read_r);
    }    

    int pair_read::read_read(pthread_mutex_t* read_mutex_1, pthread_mutex_t* read_mutex_2, pthread_mutex_t* read_mutex3, int format) {
        int status;        
        pthread_mutex_lock(read_mutex3);
        status = read1->read_read(read_mutex_1, format) * read2->read_read(read_mutex_2, format );
        pthread_mutex_unlock(read_mutex3);
        return status;
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

    int pair_read::ns_max_n(int ns_max_n) {
        int hit = read1->ns_max_n(ns_max_n) + read2->ns_max_n(ns_max_n);
        pair_read::auto_set_read_status();
        return hit;
    }

    int pair_read::min_qual_score(int min_qual) {
        int hit = read1->min_qual_score(min_qual) + read2->min_qual_score(min_qual);
        pair_read::auto_set_read_status();
        return hit;
    }

    int pair_read::min_qual_mean(int min_qual) {
        int hit = read1->min_qual_mean(min_qual) + read2->min_qual_mean(min_qual);
        pair_read::auto_set_read_status();
        return hit;
    }

    int pair_read::noiupac(void) {
        int hit = read1->noiupac() + read2->noiupac();
        pair_read::auto_set_read_status();
        return hit;
    }    


int pair_read::min_len(unsigned int len) {
    int hit = read1->min_len(len) + read2->min_len(len);
    pair_read::auto_set_read_status();
    return hit;
}    


int pair_read::max_len(unsigned int len) {
    int hit = read1->max_len(len) + read2->max_len(len);
    pair_read::auto_set_read_status();
    return hit;
}    


int pair_read::max_gc(float max_gc) {
    int hit = read1->max_gc(max_gc) + read2->max_gc(max_gc);
    pair_read::auto_set_read_status();
    return hit;
}

int pair_read::min_gc(float min_gc) {
    int hit = read1->min_gc(min_gc) + read2->min_gc(min_gc);
    pair_read::auto_set_read_status();
    return hit;
}

void pair_read::set_out_format(int format) {
    out_form=format;
}

int pair_read::entropy(float threshold) {
    int hit = read1->entropy(threshold) + read2->entropy(threshold);
    pair_read::auto_set_read_status();
    return hit;
}

int pair_read::dust(float threshold) {
    int hit = read1->dust(threshold) + read2->dust(threshold);
    pair_read::auto_set_read_status();
    return hit;
}

int pair_read::trim_tail_left(int num) {
    int hit = read1->trim_tail_left(num) + read2->trim_tail_left(num);
    pair_read::auto_set_read_status();
    return hit;
}    

int pair_read::trim_tail_right(int num) {
    int hit = read1->trim_tail_right(num) + read2->trim_tail_right(num);
    pair_read::auto_set_read_status();
    return hit;
}    

void pair_read::rm_header(void) {
    read1->rm_header();
    read2->rm_header();
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

int pair_read::trim_qual_right(string type, string rule, int step, int window_size, float threshold ) {
    int hit = read1->trim_qual_right(type, rule, step, window_size, threshold) + read2->trim_qual_right(type, rule, step, window_size, threshold);
    pair_read::auto_set_read_status();
    return hit;
} 

int pair_read::trim_qual_left(string type, string rule, int step, int window_size, float threshold ) {
    int hit = read1->trim_qual_left(type, rule, step, window_size, threshold) + read2->trim_qual_left(type, rule, step, window_size, threshold);
    pair_read::auto_set_read_status();
    return hit;
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
