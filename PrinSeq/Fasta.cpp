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
        else // TODO: need to add valid sequence checks
            IncrementBaseCount(currentLine.size());
    }
    fastaFile.close();
    PrintStats();
}

void Fasta::ProcessOptions(){
    if (optionMap.IsOptionPresent("amino")) {
        amino = optionMap.GetBoolValue("amino");
    }
    
    if (optionMap.IsOptionPresent("fasta")) {
        inputFileName = optionMap.GetStringValue("fasta") ;
        string fileType = CheckFileFormat.CheckFormat(inputFileName, amino);
        if (fileType.compare("uknown") == 0){
            cout << "Could not find input file " << '"' << inputFileName << '"'<< endl;
            return;
        }
    }
    
    if (optionMap.IsOptionPresent("seq_num")) {
        seqNum = optionMap.GetIntValue("seq_num");
    }
    
    if (optionMap.IsOptionPresent("trim_left")) {
        trimLeftAmnt = optionMap.GetIntValue("trim_left");
    }
    
    if (optionMap.IsOptionPresent("trim_right")) {
        trimRightAmnt =  optionMap.GetIntValue("trim_right");
    }
    
    if (optionMap.IsOptionPresent("trim_qual_left")) {
        trimQualLeft = optionMap.GetIntValue("trim_qual_left");
    }
    
    if (optionMap.IsOptionPresent("trim_qual_right")) {
        trimQualRight = optionMap.GetIntValue("trim_qual_right");
    }
    
    if (optionMap.IsOptionPresent("trim_tail_left")) {
        trimTailLeft = optionMap.GetIntValue("trim_tail_left");
    }
    
    if (optionMap.IsOptionPresent("trim_tail_right")) {
        trimTailRight = optionMap.GetIntValue("trim_tail_right");
    }
    
    if (optionMap.IsOptionPresent("trim_ns_left")) {
        trimNSLeft = optionMap.GetIntValue("trim_ns_left");
    }
    
    if (optionMap.IsOptionPresent("trim_ns_right")) {
        trimNSRight = optionMap.GetIntValue("trim_ns_right");
    }
    
    if (optionMap.IsOptionPresent("trim_to_len")) {
        trimToLen = optionMap.GetIntValue("trim_to_len");
    }
    
    if (optionMap.IsOptionPresent("min_len")) {
        minLength = optionMap.GetIntValue("min_len");
    }
    
    if (optionMap.IsOptionPresent("max_len")) {
        maxLength = optionMap.GetIntValue("max_len");
    }
    
    if (optionMap.IsOptionPresent("out_good")){
        goodFileName = optionMap.GetStringValue("out_good");
    }
    
    if (optionMap.IsOptionPresent("out_bad")){
        badFileName = optionMap.GetStringValue("out_bad");
    }
    
    if (optionMap.IsOptionPresent("out_format")) {
        outFormat = optionMap.GetIntValue("out_format");
    }
    
}

void Fasta::TrimSequence(){
    if (optionMap.IsOptionPresent("trim_left")) {
        FastaSeq.TrimSeqLeft(trimLeftAmnt);
    }
    
    if (optionMap.IsOptionPresent("trim_right")) {
        FastaSeq.TrimSeqRight(trimRightAmnt);
    }
    
    if (optionMap.IsOptionPresent("trim_qual_left")) {
        // TODO
    }
    
    if (optionMap.IsOptionPresent("trim_qual_right")) {
        // TODO
    }
    
    if (optionMap.IsOptionPresent("trim_tail_left")) {
        
        if(FastaSeq.GetDNASeq().find(CreateTail('A', trimTailLeft))){
            TrimTailLeft();
        }
        else if(FastaSeq.GetDNASeq().find(CreateTail('T', trimTailLeft))){
            TrimTailLeft();
        }
    }
    
    if (optionMap.IsOptionPresent("trim_tail_right")) {
        
        if(FastaSeq.GetDNASeq().find(CreateTail('A', trimTailRight))){
            TrimTailRight();
        }
        else if(FastaSeq.GetDNASeq().find(CreateTail('T', trimTailRight))){
            TrimTailRight();
        }
    }
    
    if (optionMap.IsOptionPresent("trim_ns_left")) {
        // TODO
    }
    
    if (optionMap.IsOptionPresent("trim_ns_right")) {
        // TODO
    }
    
    if (optionMap.IsOptionPresent("trim_to_len")) {
        // TODO
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
    
    if (optionMap.IsOptionPresent("stats_info")) {
        PrintStatsInfo();
    }
    
    if (optionMap.IsOptionPresent("stats_all")) {
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

