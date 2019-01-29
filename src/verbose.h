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



class verbose {
    public:
    verbose(int k,int verb);
    verbose(void);
    
    void accumulate(void);
    void print(void);
    
    int total_min_len;
    int total_max_len;
    int total_min_cg;
    int total_max_cg;
    int total_min_qual_score;
    int total_min_qual_mean;
    int total_ns_max_n;
    int total_noiupac;
    int total_derep;
    int total_lc_entropy;
    int total_lc_dust;
    int total_trim_tail_left;
    int total_trim_tail_right;
    int total_trim_qual_left;
    int total_trim_qual_right;
    int total_trim_right;
    int total_trim_left;
    int total_pair_derep;
    
    std::vector<int>* min_len;
    std::vector<int>* max_len;
    std::vector<int>* min_cg;
    std::vector<int>* max_cg;
    std::vector<int>* min_qual_score;
    std::vector<int>* min_qual_mean;
    std::vector<int>* ns_max_n;
    std::vector<int>* noiupac;
    std::vector<int>* derep;
    std::vector<int>* lc_entropy;
    std::vector<int>* lc_dust;
    std::vector<int>* trim_tail_left;
    std::vector<int>* trim_tail_right;
    std::vector<int>* trim_qual_left;
    std::vector<int>* trim_qual_right;
    std::vector<int>* trim_right;
    std::vector<int>* trim_left;
    std::vector<int>* pair_derep;
    
    int threads;
    int verbosity;
};
    
