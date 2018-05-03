#include "bloom_filter.hpp"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

#include "../config.h"

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

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>

#include "verbose.h"
verbose* verbose_vec;

pthread_mutex_t write_mutex=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t read_mutex=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t read_mutex2=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t read_mutex3=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t read_mutex4=PTHREAD_MUTEX_INITIALIZER; //derep filter

    char *forward_read_file = NULL;
    char *reverse_read_file = NULL;
    boost::filesystem::path p1,p2;
    string out_ext="fastq";
    int ii;
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
 //   int opterr = 0;
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
    int out_gz=0;
    int help=0;
    int ver=0;
    int fasta_in=0;

    std::string line;

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
        { "out_gz"          , no_argument       , &out_gz  ,  1 },
        { "threads"         , required_argument , NULL     , 22 },
        { "help"            , no_argument       , &help    ,  1 },
        { "version"         , no_argument       , &ver     ,  1 },
        { "FASTA"           , no_argument       , &fasta_in,  1 },
{0,0,0,0}
    };    
    
struct arg_struct {
    single_read * read;
    bloom_filter * filter;
    int thread_id;
};

struct arg_struct_pair {
    pair_read * read;
    bloom_filter * filter;
    int thread_id;
};

void* do_single (void * arguments);
void* do_pair (void* arguments);
void print_help(void);
void print_ver(void);

int main (int argc, char **argv)
{

    // Readin inout from the command line
    while ((c = getopt_long_only(argc, argv, "",longopts, NULL)) != -1)
        switch (c) {
            case 1:
                forward_read_file = optarg;
                p1 = optarg;
                break;
            case 2:
                reverse_read_file = optarg;
                p2=optarg;
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
            case 22: 
                threads=atoi(optarg);
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
        
    if (help) {
        print_help();
        return 0;
    }
    
    if (ver) {
        print_ver();
        return 0;
    }
    
    
/*    printf ("forward_read_file = %s ,reverse_read_file =%s\n ", forward_read_file, reverse_read_file);
    cout << "random string " << out_name << " out format " << out_format  << endl ;
    cout << "ns_max_n " << ns_max_n << endl ;
    cout  <<  "  extension()----------: " << p1.extension() << '\n';
    for (ii = optind; ii < argc; ii++)
        printf ("Non-option argument %s\n", argv[ii]);
*/


/////////// open input files

    istream * inFile_r;
    istream * inFile_f;
    ifstream *file;
    ifstream *file2;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf2;
    

    file = new ifstream(p1.native() , std::ios_base::in | std::ios_base::binary);
    file2 = new ifstream(p2.native() , std::ios_base::in | std::ios_base::binary);
    
    if (!(*file)) {
        cerr << "Error: can not open " << forward_read_file  << endl ;
        return 1;
    }
    if (reverse_read_file) {
        if (!(*file2)) {
            cerr << "Error: can not open " << reverse_read_file  << endl ;
            return 1;
        }
    } 
 
    
    
    if (p1.extension()==".gz") {
        inbuf.push(boost::iostreams::gzip_decompressor());
    }    
    
    
    if (p2.extension()==".gz") {
        inbuf2.push(boost::iostreams::gzip_decompressor());      
    } 
        

    inbuf.push(*file);
    inFile_f = new std::istream(&inbuf);
    inbuf2.push(*file2);
    inFile_r = new std::istream(&inbuf2);
    
 

//    while(std::getline(*inFile_f, line)) {
//        std::cout << line << std::endl;
//    }
//    std::cout << "bla" << line << "bla" << std::endl;
////////// open and name oupput files
    ofstream *tmp_bad_out_file_R1;
    ofstream *tmp_single_out_file_R1;
    ofstream *tmp_good_out_file_R1;
    ofstream *tmp_bad_out_file_R2;
    ofstream *tmp_single_out_file_R2;
    ofstream *tmp_good_out_file_R2;
    
    ostream *bad_out_file_R1= NULL;
    ostream *single_out_file_R1= NULL;
    ostream *good_out_file_R1= NULL;
    ostream *bad_out_file_R2= NULL;
    ostream *single_out_file_R2= NULL;
    ostream *good_out_file_R2= NULL;
    
    boost::iostreams::filtering_streambuf<boost::iostreams::output> out_bad_R1_buf;
    boost::iostreams::filtering_streambuf<boost::iostreams::output> out_bad_R2_buf;
    boost::iostreams::filtering_streambuf<boost::iostreams::output> out_single_R1_buf;
    boost::iostreams::filtering_streambuf<boost::iostreams::output> out_single_R2_buf;
    boost::iostreams::filtering_streambuf<boost::iostreams::output> out_good_R1_buf;
    boost::iostreams::filtering_streambuf<boost::iostreams::output> out_good_R2_buf;
    
    if (out_format == 1 ) { out_ext = "fasta";}
    if (out_gz == 1 ) { out_ext = out_ext + ".gz"; }
    
    if (reverse_read_file) {
        tmp_bad_out_file_R1= new std::ofstream(out_name  + "_bad_out_R1." + out_ext );
        tmp_single_out_file_R1= new std::ofstream(out_name  + "_single_out_R1." + out_ext  );
        tmp_good_out_file_R1= new std::ofstream(out_name  + "_good_out_R1." + out_ext);
        tmp_bad_out_file_R2= new std::ofstream(out_name  + "_bad_out_R2." + out_ext );
        tmp_single_out_file_R2= new std::ofstream(out_name  + "_single_out_R2." + out_ext);
        tmp_good_out_file_R2= new std::ofstream(out_name + "_good_out_R2."  + out_ext);
        if (out_gz) {
            out_bad_R1_buf.push(boost::iostreams::gzip_compressor());
            out_bad_R2_buf.push(boost::iostreams::gzip_compressor());
            out_single_R1_buf.push(boost::iostreams::gzip_compressor());
            out_single_R2_buf.push(boost::iostreams::gzip_compressor());
            out_good_R1_buf.push(boost::iostreams::gzip_compressor());
            out_good_R2_buf.push(boost::iostreams::gzip_compressor());
        }
        out_bad_R1_buf.push(*tmp_bad_out_file_R1);
        out_bad_R2_buf.push(*tmp_bad_out_file_R2);
        out_single_R1_buf.push(*tmp_single_out_file_R1);
        out_single_R2_buf.push(*tmp_single_out_file_R2);
        out_good_R1_buf.push(*tmp_good_out_file_R1);
        out_good_R2_buf.push(*tmp_good_out_file_R2);
        
        bad_out_file_R1= new std::ostream(&out_bad_R1_buf);
        bad_out_file_R2= new std::ostream(&out_bad_R2_buf);
        single_out_file_R1= new std::ostream(&out_single_R1_buf);
        single_out_file_R2= new std::ostream(&out_single_R2_buf);
        good_out_file_R1= new std::ostream(&out_good_R1_buf);
        good_out_file_R2= new std::ostream(&out_good_R2_buf);
        
        
    } else {
        tmp_good_out_file_R1= new std::ofstream(out_name  + "_good_out." + out_ext);
        tmp_bad_out_file_R1= new std::ofstream(out_name  + "_bad_out." + out_ext);
        if (out_gz) {
            out_good_R1_buf.push(boost::iostreams::gzip_compressor());
            out_bad_R1_buf.push(boost::iostreams::gzip_compressor());
        }
        out_good_R1_buf.push(*tmp_good_out_file_R1);
        out_bad_R1_buf.push(*tmp_bad_out_file_R1);
        
        good_out_file_R1= new std::ostream(&out_good_R1_buf);
        bad_out_file_R1= new std::ostream(&out_bad_R1_buf);
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
    verbose_vec= new verbose(threads);
    
    if (reverse_read_file) {
        ////////////////////////////////////////for pair end
        vector<pair_read*> v2(threads);
        vector<pthread_t> tthreads(threads);
        vector<arg_struct_pair> ttt2(threads);
        for (ii=0 ; ii<threads ; ii++){
       // v[ii].set_inputs(inFile_f);
            v2[ii] = new pair_read(*inFile_f,*inFile_r);
            v2[ii]-> set_outputs(*bad_out_file_R1,*single_out_file_R1,*good_out_file_R1,*bad_out_file_R2,*single_out_file_R2,*good_out_file_R2);
            v2[ii]-> set_out_format(out_format);
            ttt2[ii].read= v2[ii];
            ttt2[ii].filter= filter;
            ttt2[ii].thread_id=ii;
        }
        
        for (ii=0 ; ii<threads ; ii++){
            pthread_create(&tthreads[ii],NULL,do_pair, (void *) &ttt2[ii]);
        }
        
        for (ii=0 ; ii<threads ; ii++){
            pthread_join(tthreads[ii],NULL);   
        }
       
    /////////////////////////////////////////// for single end    
    } else {
             //////////// pthreads magic
        vector<single_read> v(threads,*inFile_f);
        vector<pthread_t> tthreads(threads);
        vector<arg_struct> ttt(threads);
       // declare structure for the thread
        
    
    
        for (ii=0 ; ii<threads ; ii++){
       // v[ii].set_inputs(inFile_f);
            v[ii].set_outputs(*bad_out_file_R1,*bad_out_file_R1,*good_out_file_R1);
            ttt[ii].read= & v[ii];
            ttt[ii].filter= filter;
            ttt[ii].thread_id=ii;
             
        }
        for (ii=0 ; ii<threads ; ii++){
            pthread_create(&tthreads[ii],NULL,do_single, (void *) &ttt[ii]);
        }    
    
        for (ii=0 ; ii<threads ; ii++){
            pthread_join(tthreads[ii],NULL);   
        }
    } 

//    inFile_f->close();
    if (reverse_read_file){ 
//        inFile_r->close(); 
        
    }  
    verbose_vec->accumulate();
    verbose_vec->print();
    return 0;
}


void* do_single (void * arguments) {
    struct arg_struct *args = (arg_struct*) arguments;
    single_read * read=args->read;
    bloom_filter* filter=args->filter;
    int id = args->thread_id;
    int derep_1;
    while( read->read_read( &read_mutex,fasta_in)) {
        if (trim_tail_left) {(*(verbose_vec->trim_tail_left))[id] += read->trim_tail_left(trim_tail_left);}
        if (trim_tail_right) {(*(verbose_vec->trim_tail_right))[id] += read->trim_tail_right(trim_tail_right);}
        if (trim_qual_right) {(*(verbose_vec->trim_qual_right))[id] += read->trim_qual_right("mean","lt",trim_qual_step,trim_qual_window,trim_qual_right_threshold);}
        if (trim_qual_left) {(*(verbose_vec->trim_qual_left))[id] += read->trim_qual_left("mean","lt",trim_qual_step,trim_qual_window,trim_qual_left_threshold);}
        if (ns_max_n > -1 ) {(*(verbose_vec->ns_max_n))[id] += read->ns_max_n(ns_max_n);}
        if (min_qual_mean)  {(*(verbose_vec->min_qual_mean))[id] += read->min_qual_mean(min_qual_mean);}
        if (min_qual_score) {(*(verbose_vec->min_qual_score))[id] += read->min_qual_score(min_qual_score);}
        if (noiupac) {(*(verbose_vec->noiupac))[id] += read->noiupac();}
        if (min_len) {(*(verbose_vec->min_len))[id] += read->min_len(min_len);}
        if (max_len) {(*(verbose_vec->max_len))[id] += read->max_len(max_len);}
        if (max_gc < 100) {(*(verbose_vec->max_cg))[id] += read->max_gc(max_gc);}
        if (min_gc > 0) {(*(verbose_vec->min_cg))[id] += read->min_gc(min_gc);}
        if (derep) {
            pthread_mutex_lock(& read_mutex4);
            derep_1=filter->contains(read->seq_seq);
            read->set_read_status(derep_1);
            (*(verbose_vec->derep))[id] += derep_1 ;
            filter->insert(read->seq_seq);
            pthread_mutex_unlock(& read_mutex4);
        }
        if (lc_entropy) {(*(verbose_vec->lc_entropy))[id] += read->entropy(entropy_threshold);}
        if (lc_dust) {(*(verbose_vec->lc_dust))[id] += read->dust(dust_threshold);}
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
    while(read->read_read(&read_mutex, &read_mutex2, &read_mutex3,fasta_in)) {
    int id = args->thread_id;
    int derep_1, derep_2;
        //read_rf.read1->trim_qual_right("mean","lt",5,10,30);
            if (trim_tail_left) {(*(verbose_vec->trim_tail_left))[id] += read->trim_tail_left(trim_tail_left);}
            if (trim_tail_right) {(*(verbose_vec->trim_tail_right))[id] += read->trim_tail_right(trim_tail_right);}
            if (trim_qual_right) {(*(verbose_vec->trim_qual_right))[id] += read->trim_qual_right("mean","lt",trim_qual_step,trim_qual_window,trim_qual_right_threshold);}
            if (trim_qual_left) {(*(verbose_vec->trim_qual_left))[id] += read->trim_qual_left("mean","lt",trim_qual_step,trim_qual_window,trim_qual_left_threshold);}
            if (ns_max_n > -1 ) {(*(verbose_vec->ns_max_n))[id] += read->ns_max_n(ns_max_n);}
            if (min_qual_mean)  {(*(verbose_vec->min_qual_mean))[id] += read->min_qual_mean(min_qual_mean);}
            if (min_qual_score) {(*(verbose_vec->min_qual_score))[id] += read->min_qual_score(min_qual_score);}
            if (noiupac) {(*(verbose_vec->noiupac))[id] += read->noiupac();}
            if (min_len) {(*(verbose_vec->min_len))[id] += read->min_len(min_len);}
            if (max_len) {(*(verbose_vec->max_len))[id] += read->max_len(max_len);}
            if (max_gc < 100) {(*(verbose_vec->max_cg))[id] += read->max_gc(max_gc);}
            if (min_gc > 0) {(*(verbose_vec->min_cg))[id] += read->min_gc(min_gc);}
            if (derep) {
                pthread_mutex_lock(& read_mutex4);
                derep_1=filter->contains(read->read1->seq_seq);
                derep_2=filter->contains(read->read2->seq_seq);
                read->set_read_status(derep_1,derep_2);
                (*(verbose_vec->derep))[id] += derep_1 + derep_2;
                filter->insert(read->read1->seq_seq);
                filter->insert(read->read2->seq_seq);
                pthread_mutex_unlock(& read_mutex4);
            }
        
            if (lc_entropy) {(*(verbose_vec->lc_entropy))[id] += read->entropy(entropy_threshold);}
            if (lc_dust) {(*(verbose_vec->lc_dust))[id] += read->dust(dust_threshold);}
            if (rm_header) {read->rm_header();}
            pthread_mutex_lock(& write_mutex);
            read->print();
            pthread_mutex_unlock(& write_mutex);
        } 
    pthread_exit(NULL);
}

void print_help(void) {
    std::string version = PACKAGE_VERSION;
    std::string help_msn = R"*(
        PRINSEQ++ )*" + version + R"*(
            
PRINSEQ++ is a C++ implementation of the prinseq-lite.pl program. It can be used 
to filter, reformat or trim genomic and metagenomic sequence data. It is 5X faster 
than prinseq-lite.pl and uses less RAM thanks to the use of multi-threading 
and the cboost libraries. It can read and write compressed (gzip) files, drastically 
reducing the use of hard drive.

        
Option:
    -h | -help
        Print the help page; ignore other arguments.
        
    -v | -version
        Print version; ignore other arguments.
        
    -threads 
        Nuber of threads to use. Note that if more than one thread is used, output
        sequences might not be in the same order as input sequences. (Default=1)
    
    ***** INPUT OPTIONS *****
    
    -fastq <file>
        Input file in FASTQ format. Can also read a compressed (.gz) file.
        
    -fastq2 <file>
        Input file in FASTQ format for pair-end reads. Can also read a 
        compressed (.gz) file.
        
    -FASTA
        Input is in fasta format (no quality). Note that the output format is 
        still fastq by default. Quality will be treated as 31 (A) for all bases.

    ***** OUTPUT OPTION *****
    
    -out_format <int>
        Set output format. 0 FASTQ, 1 FASTA. (Default=0)
        
    -out_name <string>
        For pair-end sequences, the output files are <string>_good_out_R1 and
        <string>_good_out_R2 for pairs where both reads pass quality control,
        <string>_single_out_R1 and <string>_single_out_R2 for read that passed
        quality control but mate didn't. <string>_bad_out_R1 and <string>_bad_out_R2  
        for reads that fail quality controls. [Default = random size 6 string] 
    
    -rm_header
        Remove the header in the 3rd line of the fastq (+header -> +). This does
        not change the header in the 1st line (@header).
        
    -out_gz 
        Write the output to a compressed file (WARNING this can be really SLOW)
        
    ***** FILTER OPTION ******
        
    -min_len <int>
        Filter sequence shorter than min_len.
    
    -max_len <int>
        Filter sequence longer than max_len.
        
    -min_gc <float>
        Filter sequence with GC percent content below min_gc.
    
    -max_gc <float>
        Filter sequence with GC percent content above min_gc.
    
    -min_qual_score <int>
        Filter sequence with at least one base with quality score below 
        min_qual_score.
        
    -min_qual_mean <int>
        Filter sequence with quality score mean below min_qual_mean.
        
    -ns_max_n <int>
        Filter sequence with more than ns_max_n Ns.
   
    -noiupac         
        Filter sequence with characters other than A, C, G, T or N.

    -derep
        Filter duplicated sequences. This only remove exact duplicates.
        
    -lc_entropy=[float]
        Filter sequences with entropy lower than [float]. [float] should be in
        the 0-1 interval. (Default=0.5)

    -lc_dust=[float]
        Filter sequences with dust_score lower than [float]. [float] should be in
        the 0-1 interval. (Default=0.5)
        
    ***** TRIM OPTIONS *****

    -trim_tail_left <integer>
        Trim poly-A/T tail with a minimum length of <integer> at the
        5'-end.

    -trim_tail_right <integer>
        Trim poly-A/T tail with a minimum length of <integer> at the
        3'-end.

    -trim_qual_rule <string>
        Rule to use to compare quality score to calculated value. Allowed
        options are lt (less than), gt (greater than) and et (equal to).
        [default: lt]

    -trim_qual_left=[float]
        Trim recursively from the 3'-end chunks of length -trim_qual_step if the
        mean quality of the first -trim_qual_window bases is less than [float]. 
        (Default=20)
        
    -trim_qual_right=[float]
        Trim recursively from the 5'-end chunks of length -trim_qual_step if the
        mean quality of the last -trim_qual_window bases is less than [float]. 
        (Default=20)    

    -trim_qual_window [int]
        Size of the window used by trim_qual_left and trim_qual_right (Default=5)

    -trim_qual_step [int]
        Step size used by trim_qual_left and trim_qual_right (Default=2)
    
    -trim_qual_type <string>
        Type of quality score calculation to use. Allowed options are min,
        mean, max and sum. [default= min]
        )*";
    std::cout << help_msn << std::endl;
    
}

void print_ver(void) {
    std::string version = PACKAGE_VERSION;
    std::string help_msn = R"*(PRINSEQ++ )*" + version;
    std::cout << help_msn << std::endl;
}
