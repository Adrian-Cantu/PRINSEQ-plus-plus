

#ifndef MATH
#define MATH
#include <math.h>
#endif

#ifndef NUMERIC
#define NUMERIC
#include <numeric>
#endif

//#include "reads.h"
//#include <unordered_map>
//#include <algorithm> 
//#include <numeric>
#include "verbose.h"

verbose::verbose(int k, int verb) : threads(k), verbosity(verb)  {
    min_len=new std::vector<int>(threads,0);
    max_len=new std::vector<int>(threads,0);
    min_cg=new std::vector<int>(threads,0);
    max_cg=new std::vector<int>(threads,0);
    min_qual_score = new std::vector<int>(threads,0);
    min_qual_mean= new std::vector<int>(threads,0);
    ns_max_n= new std::vector<int>(threads,0);
    noiupac= new std::vector<int>(threads,0);
    derep= new std::vector<int>(threads,0);
    lc_entropy= new std::vector<int>(threads,0);
    lc_dust= new std::vector<int>(threads,0);
    trim_tail_left= new std::vector<int>(threads,0);
    trim_tail_right= new std::vector<int>(threads,0);
    trim_qual_left= new std::vector<int>(threads,0);
    trim_qual_right= new std::vector<int>(threads,0);
    trim_left= new std::vector<int>(threads,0);
    trim_right= new std::vector<int>(threads,0);
    pair_derep= new std::verctor<int>(threads,0);
}

void verbose::accumulate(void) {
    total_min_len=std::accumulate((*min_len).begin(), (*min_len).end(), 0);
    total_max_len=std::accumulate((*max_len).begin(), (*max_len).end(), 0);
    total_min_cg=std::accumulate((*min_cg).begin(), (*min_cg).end(), 0);
    total_max_cg=std::accumulate((*max_cg).begin(), (*max_cg).end(), 0);
    total_min_qual_score=std::accumulate((*min_qual_score).begin(), (*min_qual_score).end(), 0);
    total_min_qual_mean=std::accumulate((*min_qual_mean).begin(),(*min_qual_mean).end(),0);
    total_ns_max_n=std::accumulate((*ns_max_n).begin(), (*ns_max_n).end(), 0);
    total_noiupac=std::accumulate((*noiupac).begin(), (*noiupac).end(), 0);
    total_derep=std::accumulate((*derep).begin(), (*derep).end(), 0);
    total_lc_entropy=std::accumulate((*lc_entropy).begin(),(*lc_entropy).end(),0);
    
    total_lc_dust=std::accumulate((*lc_dust).begin(), (*lc_dust).end(), 0);
    total_trim_tail_left=std::accumulate((*trim_tail_left).begin(), (*trim_tail_left).end(), 0);
    total_trim_tail_right=std::accumulate((*trim_tail_right).begin(), (*trim_tail_right).end(), 0);
    total_trim_qual_left=std::accumulate((*trim_qual_left).begin(), (*trim_qual_left).end(), 0);
    total_trim_qual_right=std::accumulate((*trim_qual_right).begin(), (*trim_qual_right).end(), 0);
    total_trim_left=std::accumulate((*trim_left).begin(), (*trim_left).end(), 0);
    total_trim_right=std::accumulate((*trim_right).begin(), (*trim_right).end(), 0);    
    total_pair_derep=std::accumulate((*pair_derep).begin(),(*pair_derep).end(),0)
}

void verbose::print(void){
    if (verbosity == 1 ) {
        if (total_min_len) { std::cout << total_min_len <<" reads removed by -min_len" << std::endl;}
        if (total_max_len) { std::cout << total_max_len <<" reads removed by -max_len" << std::endl;}
        if (total_min_cg) { std::cout << total_min_cg <<" reads removed by -min_cg" << std::endl;}
        if (total_max_cg) { std::cout << total_max_cg <<" reads removed by -max_cg" << std::endl;}
        if (total_min_qual_score) { std::cout << total_min_qual_score <<" reads removed by -min_qual_score" << std::endl;}
        if (total_min_qual_mean) { std::cout << total_min_qual_mean <<" reads removed by -min_qual_mean" << std::endl;}
        if (total_ns_max_n) { std::cout << total_ns_max_n <<" reads removed by -ns_max_n" << std::endl;}
        if (total_noiupac) { std::cout << total_noiupac <<" reads removed by -noiupac" << std::endl;}
        if (total_derep) { std::cout << total_derep <<" reads removed by -derep" << std::endl;}
        if (total_lc_entropy) {  std::cout << total_lc_entropy <<" reads removed by -lc_entropy" << std::endl;}
        if (total_lc_dust) { std::cout << total_lc_dust <<" reads removed by -lc_dust" << std::endl;}
        if (total_trim_tail_left) { std::cout << total_trim_tail_left <<" reads removed by -trim_tail_left" << std::endl;}
        if (total_trim_tail_right) { std::cout << total_trim_tail_right <<" reads removed by -trim_tail_right" << std::endl;}
        if (total_trim_qual_left) { std::cout << total_trim_qual_left <<" reads removed by -trim_qual_left" << std::endl;}
        if (total_trim_qual_right) { std::cout << total_trim_qual_right <<" reads removed by -trim_qual_right" << std::endl;}
        if (total_trim_left) { std::cout << total_trim_left <<" reads removed by -trim_left" << std::endl;}
        if (total_trim_right) { std::cout << total_trim_right <<" reads removed by -trim_right" << std::endl;}
        if (total_pair_derep) { std::cout << total_pair_derep <<" reads removed by -pair_derep" << std::endl;}
    } else if (verbosity==2) {
        std::cout << total_min_len        << std::endl;
        std::cout << total_max_len        << std::endl;
        std::cout << total_min_cg         << std::endl;
        std::cout << total_max_cg         << std::endl;
        std::cout << total_min_qual_score << std::endl;
        std::cout << total_min_qual_mean  << std::endl;
        std::cout << total_ns_max_n       << std::endl;
        std::cout << total_noiupac        << std::endl;
        std::cout << total_derep          << std::endl;
        std::cout << total_lc_entropy     << std::endl;
        std::cout << total_lc_dust        << std::endl;
        std::cout << total_trim_tail_left << std::endl;
        std::cout << total_trim_tail_right<< std::endl;
        std::cout << total_trim_qual_left << std::endl;
        std::cout << total_trim_qual_right<< std::endl;
        std::cout << total_trim_left << std::endl;
        std::cout << total_trim_right<< std::endl;
        std::cout << total_pair_derep << std::endl;
    }
}
