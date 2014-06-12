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
    baseCount = 0;
    outFormat = 0;
    argc = 0;
    seqNum = 0;
    trimLeft = 0;
    trimRight = 0;
    trimQualLeft = 0;
    trimQualRight = 0;
    trimTailLeft = 0;
    trimTailRight = 0;
    trimNSLeft = 0;
    trimNSRight = 0;
    trimToLen = 0;
    minLen = 0;
    maxLen = 0;
    srand((unsigned)time(0));
}

Fasta::Fasta(int count, char *fileLoc[]){
    aa = 0;  //amino acid flag
    //base = &numbase;
    seqCount = 0;
    outFormat = 1;
    fileName = "none";
    argc = count;
    string name = fileLoc[2]; // Fix later
    bOutFN = name + "_prinseq_bad_" + RandFN(); // Fix later
    gOutFN = name + "_prinseq_good_" + RandFN(); // Fix later
    
    ParseOptions(count, fileLoc);
    srand((unsigned)time(0));
    ProcessFile(fileLoc, aa);
}

void Fasta::ProcessFile(char *fileLoc[], bool aa){
    string name = *fileLoc;
    string fileType = Cff.CheckFormat(name, aa);
    string currentLine;
    string label;
    string sequence;
    if (fileType.compare("uknown") == 0){
        cout << "Could not find input file " << '"' << name << '"'<< endl;
        return;
    }
    indata.open(name);
    OutGFile.open(gOutFN);
    OutBFile.open(bOutFN);
    
    while (getline(indata, currentLine)) {
        if (currentLine[0] == '>') {
            label = currentLine;
            currentLine.clear();
        }
        else
            if (currentLine.length() < minLen) {
                /// Need to write
                currentLine.clear();
            }
    }
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
        
        /***** FILTER OPTIONS *****/
        ("min_len", po::value<int>(), "TBA")
        // Filter sequence shorter than min_len.    ("max_len", po::value<int>(), "TBA")
        // Filter sequence longer than max_len.
        
        ("range_len", po::value<string>(), "TBA")
        // Filter sequence by length range. Multiple range values should be
        // separated by comma without spaces.
        
        // Example: -range_len 50-100,250-300
        
        ("min_gc", po::value<int>(), "TBA")
        // Filter sequence with GC content below min_gc.
        
        ("max_gc", po::value<int>(), "TBA")
        // Filter sequence with GC content above max_gc.
        
        ("range_gc", po::value<string>(), "TBA")
        // Filter sequence by GC content range. Multiple range values
        // should be separated by comma without spaces.
        
        // Example: -range_gc 50-60,75-90
        
        ("min_qual_score", po::value<int>(), "TBA")
        // Filter sequence with at least one quality score below
        // min_qual_score.
        
        ("max_qual_score", po::value<int>(), "TBA")
        // Filter sequence with at least one quality score above
        // max_qual_score.
        
        ("min_qual_mean", po::value<int>(), "TBA")
        // Filter sequence with quality score mean below min_qual_mean.
        
        ("max_qual_mean", po::value<int>(), "TBA")
        // Filter sequence with quality score mean above max_qual_mean.
        
        ("ns_max_p", po::value<int>(), "TBA")
        // Filter sequence with more than ns_max_p percentage of Ns.
        
        ("ns_max_n", po::value<int>(), "TBA")
        // Filter sequence with more than ns_max_n Ns.
        
        ("noniupac", "TBA")
        // Filter sequence with characters other than A, C, G, T or N.
        
        ("seq_num", po::value<int>(), "TBA")
        // Only keep the first seq_num number of sequences (that pass all
        // other filters).
        
        ("derep", po::value<int>(), "TBA")
        // Type of duplicates to filter. Allowed values are 1, 2, 3, 4 and
        // 5. Use integers for multiple selections (e.g. 124 to use type 1,
        // 2 and 4). The order does not matter. Option 2 and 3 will set 1
        // and option 5 will set 4 as these are subsets of the other
        // option.
        
        // 1 (exact duplicate), 2 (5' duplicate), 3 (3' duplicate), 4
        // (reverse complement exact duplicate), 5 (reverse complement
        // 5'/3' duplicate)
        
        ("derep_min", po::value<int>(), "TBA")
        // This option specifies the number of allowed duplicates. If you
        // want to remove sequence duplicates that occur more than x times,
        // then you would specify x+1 as the -derep_min values. For
        // examples, to remove sequences that occur more than 5 times, you
        // would specify -derep_min 6. This option can only be used in
        // combination with -derep 1 and/or 4 (forward and/or reverse exact
        // duplicates). [default : 2]
        
        ("lc_method", po::value<string>(), "TBA")
        // Method to filter low complexity sequences. The current options
        // are "dust" and "entropy". Use "-lc_method dust" to calculate the
        // complexity using the dust method.
        
        ("lc_threshold", po::value<int>(), "TBA")
        // The threshold value (between 0 and 100) used to filter sequences
        // by sequence complexity. The dust method uses this as maximum
        // allowed score and the entropy method as minimum allowed value.
        
        ("custom_params", po::value<string>(), "TBA")
        // Can be used to specify additional filters. The current set of
        // possible rules is limited and has to follow the specifications
        // below. The custom parameters have to be specified within quotes
        // (either ' or ").
        
        // Please separate parameter values with a space and separate new
        // parameter sets with semicolon (;). Parameters are defined by two
        // values: (1) the pattern (any combination of the letters
        // "ACGTN"), (2) the number of repeats or percentage of occurence
        // Percentage values are defined by a number followed by the %-sign
        // (without space). If no %-sign is given, it is assumed that the
        // given number specifies the number of repeats of the pattern.
        
        // Examples: "AAT 10" (filters out sequences containing
        // AATAATAATAATAATAATAATAATAATAAT anywhere in the sequence), "T
        // 70%" (filters out sequences with more than 70% Ts in the
        // sequence), "A 15" (filters out sequences containing
        // AAAAAAAAAAAAAAA anywhere in the sequence), "AAT 10;T 70%;A 15"
        // (apply all three filters)
        
        /***** TRIM OPTIONS *****/
        ("trim_to_len", po::value<int>(), "TBA")
        // Trim all sequence from the 3'-end to result in sequence with
        // this length.
        
        ("trim_left", po::value<int>(), "TBA")
        // Trim sequence at the 5'-end by trim_left positions.
        
        ("trim_right", po::value<int>(), "TBA")
        // Trim sequence at the 3'-end by trim_right positions.
        
        ("trim_tail_left", po::value<int>(), "TBA")
        // Trim poly-A/T tail with a minimum length of trim_tail_left at
        // the 5'-end.
        
        ("trim_tail_right", po::value<int>(), "TBA")
        // Trim poly-A/T tail with a minimum length of trim_tail_right at
        // the 3'-end.
        
        ("trim_ns_left", po::value<int>(), "TBA")
        // Trim poly-N tail with a minimum length of trim_ns_left at the
        // 5'-end.
        
        ("trim_ns_right", po::value<int>(), "TBA")
        // Trim poly-N tail with a minimum length of trim_ns_right at the
        // 3'-end.
        
        ("trim_qual_left", po::value<int>(), "TBA")
        // Trim sequence by quality score from the 5'-end with this
        // threshold score.
        
        ("trim_qual_right", po::value<int>(), "TBA")
        // Trim sequence by quality score from the 3'-end with this
        // threshold score.
        
        ("trim_qual_type", po::value<string>(), "TBA")
        // Type of quality score calculation to use. Allowed options are
        // min, mean, max and sum. [default: min]
        
        ("trim_qual_rule", po::value<string>(), "TBA")
        // Rule to use to compare quality score to calculated value.
        // Allowed options are lt (less than), gt (greater than) and et
        // (equal to). [default: lt]
        
        ("trim_qual_window", po::value<int>(), "TBA")
        // The sliding window size used to calculate quality score by type.
        // To stop at the first base that fails the rule defined, use a
        // window size of 1. [default: 1]
        
        ("trim_qual_step", po::value<int>(), "TBA")
        // Step size used to move the sliding window. To move the window
        // over all quality scores without missing any, the step size
        // should be less or equal to the window size. [default: 1]
        
        /***** REFORMAT OPTIONS *****/
        ("seq_case", po::value<string>(), "TBA")
        // Changes sequence character case to upper or lower case. Allowed
        // options are "upper" and "lower". Use this option to remove
        // soft-masking from your sequences.
        
        ("dna_rna", po::value<string>(), "TBA")
        // Convert sequence between DNA and RNA. Allowed options are "dna"
        // (convert from RNA to DNA) and "rna" (convert from DNA to RNA).
        
        ("line_width", po::value<int>(), "TBA")
        // Sequence characters per line. Use 0 if you want each sequence in
        // a single line. Use 80 for line breaks every 80 characters. Note
        // that this option only applies to FASTA output files, since FASTQ
        // files store sequences without additional line breaks. [default:
        // 60]
        
        ("rm_header", "TBA")
        // Remove the sequence header. This includes everything after the
        // sequence identifier (which is kept unchanged).
        
        ("seq_id", po::value<string>(), "TBA")
        // Rename the sequence identifier. A counter is added to each
        // identifier to assure its uniqueness.
        
        // Example: "mySeq_10" will generate the IDs (in FASTA format)
        // >()mySeq_101, >()mySeq_102, >()mySeq_103, ...
        
        /***** SUMMARY STATISTIC OPTIONS *****/
        // The summary statistic values are written to STDOUT in the form:
        // "parameter_name statistic_name value" (without the quotes). For
        // example, "stats_info reads 10000" or "stats_len max 500". Only
        // one statistic is written per line and values are separated by
        // tabs.
        
        // If you specify any statistic option, no other ouput will be
        // generated. To preprocess data, do not specify a statistics
        // option.
        
        ("stats_info", "TBA")
        // Outputs basic information such as number of reads (reads) and
        // total bases (bases).
        
        ("stats_len", "TBA")
        // Outputs minimum (min), maximum (max), range (range), mean
        // (mean), standard deviation (stddev), mode (mode) and mode value
        // (modeval), and median (median) for read length.
        
        ("stats_dinuc", "TBA")
        // Outputs the dinucleotide odds ratio for AA/TT (aatt), AC/GT
        // (acgt), AG/CT (agct), AT (at), CA/TG (catg), CC/GG (ccgg), CG
        // (cg), GA/TC (gatc), GC (gc) and TA (ta).
        
        ("stats_tag", "TBA")
        // Outputs the probability of a tag sequence at the 5'-end (prob5)
        // and 3'-end (prob3) in percentage (0..100). Provides the number
        // of predefined MIDs (midnum) and the MID sequences (midseq,
        // separated by comma, only provided if midnum >() 0) that occur in
        // more than 34/100 (approx. 3%) of the reads.
        
        ("stats_dupl", "TBA")
        // Outputs the number of exact duplicates (exact), 5' duplicates
        // (5), 3' duplicates (3), exact duplicates with reverse
        // complements (exactrevcom) and 5'/3' duplicates with reverse
        // complements (revcomp), and total number of duplicates (total).
        // The maximum number of duplicates is given under the value name
        // with an additional "maxd" (e.g. exactmaxd or 5maxd).
        
        ("stats_ns", "TBA")
        // Outputs the number of reads with ambiguous base N (seqswithn),
        // the maximum number of Ns per read (maxn) and the maximum
        // percentage of Ns per read (maxp). The maxn and maxp value are
        // not necessary from the same sequence.
        
        ("stats_all", "Outputs all available summary statistics.")
        // Outputs all available summary statistics.
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
        
        if (vm.count("seq_num")) {
            seqNum = vm["seq_num"].as<int>();
        }
        
        if (vm.count("trim_left")) {
            trimLeft = vm["trim_left"].as<int>();
        }
        
        if (vm.count("trim_right")) {
            trimRight = vm["trim_right"].as<int>();
        }
        
        if (vm.count("trim_qual_left")) {
            trimQualLeft = vm["trim_qual_left"].as<int>();
        }
        
        if (vm.count("trim_qual_right")) {
            trimQualRight = vm["trim_qual_right"].as<int>();
        }
        
        if (vm.count("trim_ns_left")) {
            trimNSLeft = vm["trim_ns_left"].as<int>();
        }
        
        if (vm.count("trim_ns_right")) {
            trimNSRight = vm["trim_ns_right"].as<int>();
        }
        
        if (vm.count("trim_to_len")) {
            trimToLen = vm["trim_to_len"].as<int>();
        }
        
        if (vm.count("min_len")) {
            minLen = vm["min_len"].as<int>();
        }
        
        if (vm.count("max_len")) {
            maxLen = vm["max_len"].as<int>();
        }
        
        //////////////////////////////////////////////////
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

        if (vm.count("out_format")) {
            outFormat = vm["out_format"].as<int>();
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
   
	//aa = amino;
	
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

void Fasta::ProcessData(){
    good = 1;
    int begin = 0;
    int end = 0;
    
}



