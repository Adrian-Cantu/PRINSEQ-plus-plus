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
    SetDefaultValues();
}

Fasta::Fasta(int optionCount, char *OptionsArray[]){
    SetDefaultValues();
    string name = OptionsArray[2]; // Retrive name from command line arguemnts
    
    badFileName = name + "_prinseq_bad_" + RandFN() + ".fasta"; // Fix later
    goodFileName = name + "_prinseq_good_" + RandFN() + ".fasta"; // Fix later
    
    DefineOptions(optionCount, OptionsArray);
}

void Fasta::DefineOptions(int numberOfOptions, char *OptionsArray[]){
    //////////////////////////////////////
    // Cmd Line Descriptions            //
    //////////////////////////////////////
    try {
        // Defining and adding all allowable options 
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
        
        ("amino", po::value<bool>(), "set amino acid")
        // Input is amino acid (protein) sequences instead of nucleic acid (DNA or RNA) sequences. Allowed amino
        // acid characters: ABCDEFGHIKLMNOPQRSTUVWYZXabcdefghiklmmopqrstuvwyzx*- and allowed nucleic acid
        // characters: ACGTURYKMSWBDHVNXacgturykmswbdhvnx-
        //
        // The following options are ignored for -amino: stats_dinuc,stats_tag,stats_ns,dna_rna
        
        
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
        // Filter sequence shorter than min_len.
        
        ("max_len", po::value<int>(), "TBA")
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
        
        // Examples: "aminoT 10" (filters out sequences containing
        // aminoTaminoTaminoTaminoTaminoTaminoTaminoTaminoTaminoTaminoT anywhere in the sequence), "T
        // 70%" (filters out sequences with more than 70% Ts in the
        // sequence), "A 15" (filters out sequences containing
        // aminoaminoaminoaminoaminoaminoaminoA anywhere in the sequence), "aminoT 10;T 70%;A 15"
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
        // Outputs the dinucleotide odds ratio for amino/TT (aminott), AC/GT
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
        
        //        po::variables_map vm; // Holds all options from cmd line
        po::store(po::parse_command_line(numberOfOptions, OptionsArray, desc), vm); // stores options in vm
        po::notify(vm);
        
    }
    catch(std::exception& e) {
        cerr << "error: " << e.what() << "\n";
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }
}

void Fasta::ProcessFile(){
    ProcessOptions();
    string currentLine;
    
    fastaFile.open(inputFileName);
    while (getline(fastaFile, currentLine)) {

        trim(currentLine);
        if (currentLine[0] == '>') {
            FastaSeq.SetID(currentLine);
        }
        else{
            FastaSeq.SetDNA(currentLine);
            ApplyFilters();
        }
        if (IsSeqID(currentLine)) {
            IncrementSeqCount();
        }
        else // need to add valid sequence checks
            IncrementBaseCount(currentLine.size());
    }
    fastaFile.close();
    PrintStats();
}

void Fasta::ProcessOptions(){
    
    ///////////////////////////////////////
    // Functions for cmd line parameters //
    ///////////////////////////////////////
    
    if (vm.count("amino")) {
        amino = vm["amino"].as<bool>();
    }
    
    if (vm.count("fasta")) {
        inputFileName = vm["fasta"].as<string>();
        string fileType = CheckFileFormat.CheckFormat(inputFileName, amino);
        if (fileType.compare("uknown") == 0){
            cout << "Could not find input file " << '"' << vm["fasta"].as<string>() << '"'<< endl;
            return;
        }
    } // Calls File Format Class: FormatCheck.cpp
    
    if (vm.count("seq_num")) {
        seqNum = vm["seq_num"].as<int>();
    }
    
    if (vm.count("trim_left")) {
        trimLeftAmnt = vm["trim_left"].as<int>();
    }
    
    if (vm.count("trim_right")) {
        trimRightAmnt =  vm["trim_right"].as<int>();
    }
    
    if (vm.count("trim_qual_left")) {
        trimQualLeft = vm["trim_qual_left"].as<int>();
    }
    
    if (vm.count("trim_qual_right")) {
        trimQualRight = vm["trim_qual_right"].as<int>();
    }
    
    if (vm.count("trim_tail_left")) {
        trimTailLeft = vm["trim_tail_left"].as<int>();
    }
    
    if (vm.count("trim_tail_right")) {
        trimTailRight = vm["trim_tail_right"].as<int>();
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
        minLength = vm["min_len"].as<int>();
    }
    
    if (vm.count("max_len")) {
        maxLength = vm["max_len"].as<int>();
    }
    
    //////////////////////////////////////////////////
    if (vm.count("out_good")){
        goodFileName = vm["out_good"].as<string>();
    }
    
    if (vm.count("out_bad")){
        badFileName = vm["out_bad"].as<string>();
    }
    
    if (vm.count("out_format")) {
        outFormat = vm["out_format"].as<int>();
    }
    
}

void Fasta::TrimSequence(){
    if (vm.count("trim_left")) {
        FastaSeq.TrimSeqLeft(trimLeftAmnt);
    }
    
    if (vm.count("trim_right")) {
        FastaSeq.TrimSeqRight(trimRightAmnt);
    }
    
    if (vm.count("trim_qual_left")) {
        
    }
    
    if (vm.count("trim_qual_right")) {
        
    }
    
    if (vm.count("trim_ns_left")) {
        
    }
    
    if (vm.count("trim_ns_right")) {
        
    }
    
    if (vm.count("trim_to_len")) {
        
    }

}

void Fasta::TrimQualLeft(){
    // Need to Code Qual First
}

void Fasta::TrimQualRight(){
    // Need to Code Qual First
}

void Fasta::TrimTailLeft(){
    int trimValue = 0;
    while (FastaSeq.GetDNASeq()[trimValue] != 'A' || 'T' || 'N') {
        trimValue++;
    }
    FastaSeq.TrimSeqLeft(trimValue);
}

void Fasta::TrimTailRight(){
    long trimValue = FastaSeq.GetSeqLength();
    while (FastaSeq.GetDNASeq()[trimValue] != 'A' || 'T' || 'N') {
        trimValue--;
    }
    FastaSeq.TrimSeqRight((int)trimValue);
}

void Fasta::TrimNSLeft(){
    int trimValue = 0;
    while (FastaSeq.GetDNASeq()[trimValue] != 'N') {
        trimValue++;
    }
    FastaSeq.TrimSeqLeft(trimValue);
}

void Fasta::TrimNSRight(){
    long trimValue = FastaSeq.GetSeqLength();
    while (FastaSeq.GetDNASeq()[trimValue] != 'N') {
        trimValue--;
    }
    FastaSeq.TrimSeqRight((int)trimValue);
}

void Fasta::ApplyFilters(){
    if (minLength > 0) {
        MinLengthFilter();
    }
    else if (maxLength > 0) {
        MaxLengthFilter();
    }
}

void Fasta::MinLengthFilter(){
    if(minLength <= FastaSeq.GetSeqLength()){
        IncrementGoodSeqCount();
        IncrementGoodBaseCount();
        //WriteToGood();
    }
    else{
        IncrementBadSeqCount();
        IncrementBadBaseCount();
        //WriteToBad();
    }
}

void Fasta::MaxLengthFilter(){
    if(maxLength >= FastaSeq.GetSeqLength()){
        //IncrementGoodSeqCount();
        IncrementGoodBaseCount();
        //WriteToGood();
    }
    else{
        IncrementBadSeqCount();
        IncrementBadBaseCount();
        //WriteToBad();
    }
}

string Fasta::RandFN(){
    string filename;
    filename = to_string(rand() % 10000);
    
    return filename;
}


void Fasta::SetOutputFormat(int format){
    outFormat = format;
    cout << outFormat << endl;
}

//****** INCREMENT Seq/Base counts ******//
void Fasta::IncrementSeqCount(){
    seqCount++;
}
void Fasta::IncrementBaseCount(long size){
    baseCount += size;
}


//****** INCREMENT BAD Seq/Base counts ******//
void Fasta::IncrementBadSeqCount(){
    badSeqCount++;
}

void Fasta::IncrementBadBaseCount(){
    badBaseCount += FastaSeq.GetDNASeq().size();
}


//****** INCREMENT GOOD Seq/Base counts ******//
void Fasta::IncrementGoodSeqCount(){
    goodSeqCount = goodSeqCount + 1;
}
void Fasta::IncrementGoodBaseCount(){
    goodBaseCount += FastaSeq.GetDNASeq().size();
}


//****** GET Seq/Base counts ******//
long Fasta::GetSeqCount(){
    return seqCount;
}

long Fasta::GetBaseCount(){
    return baseCount;
}


//****** GET BAD Seq/Base counts ******//
long Fasta::GetBadSeqCount(){
    return badSeqCount;
}

long Fasta::GetBadBaseCount(){
    return badBaseCount;
}


//****** GET GOOD Seq/Base counts ******//
long Fasta::GetGoodSeqCount(){
    return goodSeqCount;
}

long Fasta::GetGoodBaseCount(){
    return goodBaseCount;
}


//****** WRITE Files ******//
void Fasta::WriteToGood()
{
    GoodFileStream.open(goodFileName, ios::app);
    GoodFileStream << FastaSeq.GetID() << endl;
    GoodFileStream << FastaSeq.GetDNASeq() << endl;
    GoodFileStream.close();
}

void Fasta::WriteToBad()
{
    BadFileStream.open(badFileName, ios::app);
    BadFileStream << FastaSeq.GetID() << endl;
    BadFileStream << FastaSeq.GetDNASeq() << endl;
    BadFileStream.close();
}

//****** PRINT Stats ******//
void Fasta::PrintStats(){
    PrintStandardStats();
    
    if (vm.count("stats_info")) {
        PrintStatsInfo();
    }
    
    if (vm.count("stats_all")) {
        PrintStats_All();
    }
}

void Fasta::PrintStandardStats(){
    double mean = double(GetBaseCount())/double(GetSeqCount());
    cout << "Input and filter stats:" << endl;
    cout << "\t\tInput sequences: " << GetSeqCount() << endl;
    cout << "\t\tInput bases: " << GetBaseCount() << endl;
    cout << "\t\tInput mean length: " << fixed << setprecision(2) << showpoint << mean << endl;
    cout << "\t\tGood sequences: " << GetGoodSeqCount() << endl;
    cout << "\t\tGood bases: " << GetGoodBaseCount() << endl;
    cout << "\t\tBad sequences: " << GetBadSeqCount() << endl;
    cout << "\t\tBad bases: " << GetBadBaseCount() << endl;
    cout << "\t\tSequences filtered by specified parameters: "<< endl << "\t\tnone" << endl;
}

void Fasta::PrintStatsInfo(){
    cout << "\t\tInput sequences: " << GetSeqCount() << endl;
    cout << "\t\tInput bases: " << GetBaseCount() << endl;
}

void Fasta::PrintStats_All(){
    cout << "\t\tInput sequences: " << GetSeqCount() << endl;
    cout << "\t\tInput bases: " << GetBaseCount() << endl;
}

bool Fasta::IsSeqID(string s){
    static const regex e1("^>(\\S+)\\s*(.*)$");
    return regex_match(s,e1);
}


void Fasta::SetDefaultValues(){
    amino = 1;
    seqCount = 0;
    baseCount = 0;
    badSeqCount = 0;
    badBaseCount = 0;
    goodSeqCount = 0;
    goodBaseCount = 0;
    outFormat = 1;
    inputFileName = "none";
    outFormat = 0;
    seqNum = 0;
    trimLeftAmnt = 0;
    trimRightAmnt = 0;
    trimQualLeft = 0;
    trimQualRight = 0;
    trimTailLeft = 0;
    trimTailRight = 0;
    trimNSLeft = 0;
    trimNSRight = 0;
    trimToLen = 0;
    minLength = 0;
    maxLength = 0;
}

