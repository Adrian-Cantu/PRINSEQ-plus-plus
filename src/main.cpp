#include "bloom_filter.hpp"
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

#ifndef PTHREAD
#define PTHREAD
#include <pthread.h>
#endif


pthread_mutex_t write_mutex=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t read_mutex=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t read_mutex2=PTHREAD_MUTEX_INITIALIZER;

    char *forward_read_file = NULL;
    char *reverse_read_file = NULL;
    string out_ext="fastq";
    int index;
    int out_format=0; // 0 fastq, 1 fasta
    int min_qual_score=0;
    int ns_max_n=-1;
    int min_qual_mean=0;
    int noiupac=0;
    int c;
    int min_len=0;
    int max_len=0;
    float max_gc=100;
    float min_gc=0;
    int opterr = 0;
    int derep;
    float entropy_threshold=0.5 ;
    int lc_entropy=0;
    float dust_threshold=0.5 ;
    int lc_dust=0;
    
    float trim_qual_right_threshold=20;
    int trim_qual_right=0;
    float trim_qual_left_threshold=20;
    int trim_qual_left=0;
    string trim_qual_type="min";
    string trim_qual_rule="lt";
    int trim_qual_window=5;
    int trim_qual_step=2;
    string out_name=random_string(6);
    int rm_header=0;
    int trim_tail_left=0;
    int trim_tail_right=0;
    int threads=5;
    int ii=0;

    struct option longopts[] = {
        { "fastq"           , required_argument , NULL     ,  1 },
        { "fastq2"          , required_argument , NULL     ,  2 },
        { "out_format"      , required_argument , NULL     ,  3 },
        { "min_qual_score"  , required_argument , NULL     ,  4 },
        { "ns_max_n"        , required_argument , NULL     ,  5 },
        { "min_qual_mean"   , required_argument , NULL     ,  6 },
        { "noiupac"         , no_argument       , &noiupac ,  1 },
        { "derep"           , no_argument       , &derep   ,  1 },
        { "min_len"         , required_argument , NULL     ,  7 },
        { "max_len"         , required_argument , NULL     ,  8 },
        { "max_gc"          , required_argument , NULL     ,  9 },
        { "min_gc"          , required_argument , NULL     , 10 },
        { "lc_entropy"      , optional_argument , NULL     , 11 },
        { "lc_dust"         , optional_argument , NULL     , 12 },
        { "trim_qual_right" , optional_argument , NULL     , 13 },
        { "trim_qual_left"  , optional_argument , NULL     , 14 },
        { "trim_qual_type"  , required_argument , NULL     , 15 },
        { "trim_qual_rule"  , required_argument , NULL     , 16 },
        { "trim_qual_window", required_argument , NULL     , 17 },
        { "trim_qual_step"  , required_argument , NULL     , 18 },
        { "out_name"        , required_argument , NULL     , 19 },
        { "rm_header"       , no_argument       , &rm_header, 1 },
        { "trim_tail_left"  , required_argument , NULL     , 20 },
        { "trim_tail_right" , required_argument , NULL     , 21 },
{0,0,0,0}
    };    
    
struct arg_struct {
    single_read * read;
    bloom_filter * filter;
};

struct arg_struct_pair {
    pair_read * read;
    bloom_filter * filter;
};

void* do_single (void * arguments);
void* do_pair (void* arguments);

int main (int argc, char **argv)
{

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
            case 7:
                min_len = atoi(optarg);
                break;
            case 8:
                max_len = atoi(optarg);
                break;
            case 9:
                max_gc = atof(optarg);
                break;
            case 10:
                min_gc = atoi(optarg);
                break;
            case 11:
                if (optarg != NULL) {
                    entropy_threshold = atof(optarg);
                }
                lc_entropy=1;
                break;
            case 12:
                if (optarg != NULL) {
                    dust_threshold = atof(optarg);
                }
                lc_dust=1;
                break;
            case 13:
                if (optarg != NULL) {
                    trim_qual_right_threshold = atof(optarg);
                }
                trim_qual_right=1;
                break;
            case 14:
                if (optarg != NULL) {
                    trim_qual_left_threshold = atof(optarg);
                }
                trim_qual_left=1;
                break;
            case 15:
                trim_qual_type=optarg;
                break;
            case 16:
                trim_qual_rule=optarg;
                break;
            case 17:
                trim_qual_window=atoi(optarg);
                break;
            case 18:
                trim_qual_step=atoi(optarg);
                break;
            case 19:
                out_name=optarg;
                break;
            case 20:
                trim_tail_left=atoi(optarg);
                break;
            case 21: 
                trim_tail_right=atoi(optarg);
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

    printf ("forward_read_file = %s ,reverse_read_file =%s\n ", forward_read_file, reverse_read_file);
    cout << "random string " << out_name << " out format " << out_format  << endl ;
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
    if (reverse_read_file) {
        inFile_r.open(reverse_read_file);
        if (!inFile_r) {
            cerr << "Error: can not open " << reverse_read_file  << endl ;
            return 1;
        }
    }    
      
////////// open and name oupput files
    ofstream bad_out_file_R1;
    ofstream single_out_file_R1;
    ofstream good_out_file_R1;
    ofstream bad_out_file_R2;
    ofstream single_out_file_R2;
    ofstream good_out_file_R2;

    if (out_format == 1 ) { out_ext = "fasta";}
    
    
    if (reverse_read_file) {
        bad_out_file_R1.open(out_name  + "_bad_out_R1." + out_ext );
        single_out_file_R1.open(out_name  + "_single_out_R1." + out_ext  );
        good_out_file_R1.open(out_name  + "_good_out_R1." + out_ext);
        bad_out_file_R2.open(out_name  + "_bad_out_R2." + out_ext );
        single_out_file_R2.open(out_name  + "_single_out_R2." + out_ext);
        good_out_file_R2.open(out_name + "_good_out_R2."  + out_ext);
    } else {
        bad_out_file_R1.open(out_name  + "_bad_out." + out_ext );
        good_out_file_R1.open(out_name  + "_good_out." + out_ext);
    }
    

    single_read read_f(inFile_f);
    single_read read_r(inFile_r);
    pair_read read_rf(inFile_f,inFile_r);
    if (reverse_read_file) {
        read_rf.set_outputs(bad_out_file_R1,single_out_file_R1,good_out_file_R1,
            bad_out_file_R2,single_out_file_R2,good_out_file_R2);
        read_rf.set_out_format(out_format);
    } else {
        read_f.set_outputs(bad_out_file_R1,single_out_file_R1,good_out_file_R1);
    }
    bloom_parameters parameters;
    bloom_filter *filter=NULL;
    if (derep) {    
        parameters.projected_element_count = 10000000;
        parameters.false_positive_probability = 0.000001; // 1 in 10000
        parameters.random_seed = 0xA5A5A5A5;
        if (!parameters) {
            std::cout << "Error - Invalid set of bloom filter parameters!" << std::endl;
            return 1;
        }
        parameters.compute_optimal_parameters();
        filter  = new bloom_filter(parameters);
    }
 
 
    
 
 
 
    // main loop
    if (reverse_read_file) {
        ////////////////////////////////////////for pair end
        vector<pair_read> v2(threads);
        vector<pthread_t> tthreads(threads);
        vector<arg_struct_pair> ttt2(threads);
        for (ii=0 ; ii<threads ; ii++){
       // v[ii].set_inputs(inFile_f);
            v2[ii].set_inputs(inFile_f,inFile_r);
            v2[ii].set_outputs(bad_out_file_R1,single_out_file_R1,good_out_file_R1,bad_out_file_R2,single_out_file_R2,good_out_file_R2);
            v2[ii].set_out_format(out_format);
            ttt2[ii].read= & v2[ii];
            ttt2[ii].filter= filter;
            pthread_create(&tthreads[ii],NULL,do_pair, (void *) &ttt2[ii]); 
        }
        
        for (ii=0 ; ii<threads ; ii++){
            pthread_join(tthreads[ii],NULL);   
        }
       
    /////////////////////////////////////////// for single end    
    } else {
             //////////// pthreads magic
        vector<single_read> v(threads,inFile_f);
        vector<pthread_t> tthreads(threads);
        vector<arg_struct> ttt(threads);
       // declare structure for the thread
        
    
    
        for (ii=0 ; ii<threads ; ii++){
       // v[ii].set_inputs(inFile_f);
            v[ii].set_outputs(bad_out_file_R1,single_out_file_R1,good_out_file_R1);
            ttt[ii].read= & v[ii];
            ttt[ii].filter= filter;
            pthread_create(&tthreads[ii],NULL,do_single, (void *) &ttt[ii]); 
        }
    
        for (ii=0 ; ii<threads ; ii++){
            pthread_join(tthreads[ii],NULL);   
        }
    } 
    
    inFile_f.close();
    if (reverse_read_file){ 
        inFile_r.close(); 
        
    }  
    return 0;
}


void* do_single (void * arguments) {
    struct arg_struct *args = (arg_struct*) arguments;
    single_read * read=args->read;
    bloom_filter* filter=args->filter;
    while( read->read_read( &read_mutex)) {
        if (trim_tail_left) {read->trim_tail_left(trim_tail_left);}
        if (trim_tail_right) {read->trim_tail_right(trim_tail_right);}
        if (trim_qual_right) {read->trim_qual_right("mean","lt",trim_qual_step,trim_qual_window,trim_qual_right_threshold);}
        if (trim_qual_left) {read->trim_qual_left("mean","lt",trim_qual_step,trim_qual_window,trim_qual_left_threshold);}
        if (ns_max_n > -1 ) {read->ns_max_n(ns_max_n);}
        if (min_qual_mean)  {read->min_qual_mean(min_qual_mean);}
        if (min_qual_score) {read->min_qual_score(min_qual_score);}
        if (noiupac) {read->noiupac();}
        if (min_len) {read->min_len(min_len);}
        if (max_len) {read->max_len(max_len);}
        if (max_gc < 100) {read->max_gc(max_gc);}
        if (min_gc > 0) {read->min_gc(min_gc);}
        if (derep) {
            if(filter->contains(read->seq_seq)) { read->set_read_status(2);}
            filter->insert(read->seq_seq);
            
        }
    
        if (lc_entropy) {read->entropy(entropy_threshold);}
        if (lc_dust) {read->dust(dust_threshold);}
        if (rm_header) {read->rm_header();}
        pthread_mutex_lock(& write_mutex);
        read->print(out_format);
        pthread_mutex_unlock(& write_mutex);
    }
    pthread_exit(NULL);
    
}

void* do_pair (void * arguments) {
    struct arg_struct_pair *args = (arg_struct_pair*) arguments;
    pair_read * read=args->read;
    bloom_filter* filter=args->filter;
    while(read->read_read(&read_mutex, &read_mutex2)) {
        //read_rf.read1->trim_qual_right("mean","lt",5,10,30);
            if (trim_tail_left) {read->trim_tail_left(trim_tail_left);}
            if (trim_tail_right) {read->trim_tail_right(trim_tail_right);}
            if (trim_qual_right) {read->trim_qual_right("mean","lt",trim_qual_step,trim_qual_window,trim_qual_right_threshold);}
            if (trim_qual_left) {read->trim_qual_left("mean","lt",trim_qual_step,trim_qual_window,trim_qual_left_threshold);}
            if (ns_max_n > -1 ) {read->ns_max_n(ns_max_n);}
            if (min_qual_mean)  {read->min_qual_mean(min_qual_mean);}
            if (min_qual_score) { read->min_qual_score(min_qual_score);}
            if (noiupac) {read->noiupac();}
            if (min_len) {read->min_len(min_len);}
            if (max_len) {read->max_len(max_len);}
            if (max_gc < 100) {read->max_gc(max_gc);}
            if (min_gc > 0) {read->min_gc(min_gc);}
            if (derep) {
                read->set_read_status(filter->contains(read->read1->seq_seq),filter->contains(read->read2->seq_seq));
                filter->insert(read->read1->seq_seq);
                filter->insert(read->read2->seq_seq);
            }
        
            if (lc_entropy) {read->entropy(entropy_threshold);}
            if (lc_dust) {read->dust(dust_threshold);}
            if (rm_header) {read->rm_header();}
            pthread_mutex_lock(& write_mutex);
            read->print();
            pthread_mutex_unlock(& write_mutex);
        } 
    pthread_exit(NULL);
}    