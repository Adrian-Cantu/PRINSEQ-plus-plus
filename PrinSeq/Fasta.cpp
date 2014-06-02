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
    //base = &numbase;
    seqCount = 0;
    outFormat = 1;
    fileName = "none";
    srand((unsigned)time(0));
}

Fasta::Fasta(int count, char *fileLoc[]){
    aa = 0;  //amino acid flag
    //base = &numbase;
    seqCount = 0;
    outFormat = 0;
    fileName = "none";
    argc = count;
    string name = fileLoc[2]; // Fix later
    bOutFN = name + "_prinseq_bad_" + RandFN(); // Fix later
    gOutFN = name + "_prinseq_good_" + RandFN(); // Fix later
    
    ParseOptions(count, fileLoc);
    srand((unsigned)time(0));
    
}

void Fasta::ParseOptions(int count, char *fileLoc[]){
    //////////////////////////////////////
    // Cmd Line Descriptions            //
    //////////////////////////////////////
    try {
        
        po::options_description desc("Allowed options");
        desc.add_options()
        
        ("fasta", po::value<string>(), "file type")
        // Input file in FASTA format that contains the sequence data. Use stdin instead of a file name to read
        // from STDIN (-fastq stdin). This can be useful to process compressed files using Unix pipes.
        
        ("qual", po::value<string>(), "quality file")
        // Input file in QUAL format that contains the quality data.
        
        
        /*("fastq", po::value<string>(), "file type")
         // Input file in FASTQ format that contains the sequence and quality data. Use stdin instead of a file
         // name to read from STDIN (-fasta stdin). This can be useful to process compressed files using Unix
         // pipes. */
        
        ("aa", po::value<bool>(), "set amino acid")
        // Input is amino acid (protein) sequences instead of nucleic acid (DNA or RNA) sequences. Allowed amino
        // acid characters: ABCDEFGHIKLMNOPQRSTUVWYZXabcdefghiklmmopqrstuvwyzx*- and allowed nucleic acid
        // characters: ACGTURYKMSWBDHVNXacgturykmswbdhvnx-
        //
        // The following options are ignored for -aa: stats_dinuc,stats_tag,stats_ns,dna_rna
        
        
        ("out_good", po::value<string>(), "Good output filename")
        // By default, the output files are created in the same directory as the input file containing the
        // sequence data with an additional "_prinseq_good_XXXX" in their name (where XXXX is replaced by random
        // characters to prevent overwriting previous files). To change the output filename and location,
        // specify the filename using this option. The file extension will be added automatically. Use "
        // out_good null" to prevent the program from generating the output file(s) for data passing all
        // filters. Use "-out_good stdout" to write data passing all filters to STDOUT (only for FASTA or FASTQ
        // output files).
        //
        // Example: use "file_passed" to generate the output file file_passed.fasta in the current directory
        
        ("out_bad", po::value<string>(), "Bad output filename")
        // By default, the output files are created in the same directory as the input file containing the
        // sequence data with an additional "_prinseq_bad_XXXX" in their name (where XXXX is replaced by random
        // characters to prevent overwriting previous files). To change the output filename and location,
        // specify the filename using this option. The file extension will be added automatically. Use "-out_bad
        // null" to prevent the program from generating the output file(s) for data not passing any filter. Use
        // "-out_bad stdout" to write data not passing any filter to STDOUT (only for FASTA or FASTQ output
        // files).
        //
        // Example: use "file_filtered" to generate the output file file_filtered.fasta in the current directory
        //
        // Example: "-out_good stdout -out_bad null" will write data passing filters to STDOUT and data not
        // passing any filter will be ignored
        
        
        ("out_format", po::value<int>(), "Format of output file")
        // To change the output format, use one of the following options. If not defined, the output format will
        // be the same as the input format.
        //
        // 1 (FASTA only), 2 (FASTA and QUAL), 3 (FASTQ), 4 (FASTQ and FASTA), or 5 (FASTQ, FASTA and QUAL)
        
        ("stats_all", "GOutputs all available summary statistics.")
        ;
        
        
        
        po::variables_map vm; // Holds all options from cmd line
        po::store(po::parse_command_line(argc, fileLoc, desc), vm); // stores options in vm
        po::notify(vm);
        
        ///////////////////////////////////////
        // Functions for cmd line parameters //
        ///////////////////////////////////////
        
        if (vm.count("aa")) {
            aa = vm["aa"].as<bool>();
        }// Sets AA value if specified, otherwise 0.
        
        if (vm.count("fasta")) {
            fileName = vm["fasta"].as<string>();
            string fileType = Cff.CheckFormat(vm["fasta"].as<string>(), aa);
            if (fileType.compare("uknown") == 0){
                cout << "Could not find input file " << '"' << vm["fasta"].as<string>() << '"'<< endl;
                return;
            }
        } // Calls File Format Class: FormatCheck.cpp
        
        if (vm.count("out_format")) {
            outFormat = vm["out_format"].as<int>();
        }
        
        if (vm.count("out_good")){
            WriteGood(vm["out_good"].as<string>(), aa);
        }
            
        if (vm.count("out_bad")){
            WriteBad(vm["out_bad"].as<string>(), aa);
        }
        
        if (vm.count("stats_all")) {
            cout << "Fasta Class File Name: " << fileName << endl;
            cout << "test stats all." << endl;
            Stats_All();
        }
        

    }
    catch(std::exception& e) {
        cerr << "error: " << e.what() << "\n";
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }
}

string Fasta::RandFN(){
    string filename;
    filename = rand() % 1000;
    
    return filename;
}

void Fasta::WriteGood(string filename, bool amino){
    
    gOutFN = filename+ ".fasta";
    cout << gOutFN << endl;
}
void Fasta::WriteBad(string filename, bool amino){
    
    bOutFN = filename + ".fasta";
    cout << bOutFN << endl;
}

void Fasta::SetOutputFormat(int format){
    outFormat = format;
    cout << outFormat << endl;
}

void Fasta::AddSeqCount(){
    seqCount++;
}
void Fasta::AddBaseCount(int size){
    baseCount += size;
}
int Fasta::GetBaseCount(){
    return baseCount;
}
int Fasta::GetSeqCount(){
    return seqCount;
}
void Fasta::Stats_All(){
    indata.open(fileName);
    //cout << filename << endl;
	//aa = amino;
	//cout << "QualAlarm3" << endl;
	if(!indata) { // file couldn't be opened
		string error = "File could not be opened";
		return ;
	}
    
    string lineA;
    while (getline(indata, lineA)) {
        if (ValExp(lineA)) {
            AddSeqCount();
        }
        else // need to add valid sequence checks
            AddBaseCount(lineA.size());
    }
    cout << "Sequence Count: " << seqCount << endl;
    cout << "Base Count: " << baseCount << endl;
}

bool Fasta::ValExp(string s){
    static const regex e1("^>(\\S+)\\s*(.*)$");
    return regex_match(s,e1);
}

