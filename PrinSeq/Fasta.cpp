//
//  Fasta.cpp
//  PrinSeqCPP
//
// 
//  Copyright (c) 2014 Edwards Lab. All rights reserved.
//

#include "Fasta.h"
using namespace std;
Fasta::Fasta(){
    
    aa = 0;  //amino acid flag
    base = &numbase;
    seqCount = 0;
    
}

Fasta::Fasta(int argc, char *fileLoc[]){
    
    aa = 0;  //amino acid flag
    base = &numbase;
    seqCount = 0;
    
    
    try {
        
        po::options_description desc("Allowed options");
        desc.add_options()
        ("fasta", po::value<string>(), "file type")
        ("qual", po::value<string>(), "quality file")
        ("aa", po::value<bool>(), "set amino acid")
        ("out_good", po::value<string>(), "Good output filename")
        ("out_bad", po::value<string>(), "Bad output filename")
        ("out_format", po::value<int>(), "Format of output file")
        ;
        
        po::variables_map vm;
        po::store(po::parse_command_line(argc, fileLoc, desc), vm);
        po::notify(vm);
        
        if (vm.count("fasta")) {
            string test = Cff.CheckFormat(vm["fasta"].as<string>(), vm["aa"].as<bool>());
        }
        
        if (vm.count("aa")) {
            cout << "AA was set to "
            << vm["aa"].as<bool>() << ".\n";
        } else {
            cout << "AA was not set.\n";
        }
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }
}



