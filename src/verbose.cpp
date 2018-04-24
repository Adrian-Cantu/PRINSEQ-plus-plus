

#ifndef MATH
#define MATH
#include <math.h>
#endif


//#include "reads.h"
//#include <unordered_map>
//#include <algorithm> 
//#include <numeric>
#include "verbose.h"

verbose::verbose(int k) : threads(k) {
    min_len=new std::vector<int>(threads,0);
    max_len=new std::vector<int>(threads,0);
    min_cg=new std::vector<int>(threads,0);
    max_cg=new std::vector<int>(threads,0);
}

void verbose::accumulate(void) {
    total_min_len=std::accumulate((*min_len).begin(), (*min_len).end(), 0);
}

void verbose::print(void){
    std::cout << total_min_len <<" read removed by -min_len" << std::endl;
}