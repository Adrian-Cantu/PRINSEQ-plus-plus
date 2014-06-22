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
    optionMap.DefineOptions(numberOfOptions, OptionsArray);
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
    
    if (optionMap.IsOptionPresent("amino")) {
        amino = optionMap.GetBoolValue("amino");
    }
    
    if (optionMap.IsOptionPresent("fasta")) {
        inputFileName = vm["fasta"].as<string>();
        string fileType = CheckFileFormat.CheckFormat(inputFileName, amino);
        if (fileType.compare("uknown") == 0){
            cout << "Could not find input file " << '"' << vm["fasta"].as<string>() << '"'<< endl;
            return;
        }
    } // Calls File Format Class: FormatCheck.cpp
    
    if (vm.count("seq_num")) {
        seqNum = optionMap.GetIntValue("seq_num");
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
    
    if (vm.count("trim_tail_left")) {
        
        if(FastaSeq.GetDNASeq().find(CreateTail('A', trimTailLeft))){
            TrimTailLeft();
        }
        else if(FastaSeq.GetDNASeq().find(CreateTail('T', trimTailLeft))){
            TrimTailLeft();
        }
    }
    
    if (vm.count("trim_tail_right")) {
        
        if(FastaSeq.GetDNASeq().find(CreateTail('A', trimTailRight))){
            TrimTailRight();
        }
        else if(FastaSeq.GetDNASeq().find(CreateTail('T', trimTailRight))){
            TrimTailRight();
        }
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

string Fasta::CreateTail(char ATN, int tailLength){ //create a tail consistion of A, T, or N
    ostringstream newTail;
    for (int i = 0; i < tailLength; i++) {
        newTail << ATN;
    }
    return newTail.str();
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

