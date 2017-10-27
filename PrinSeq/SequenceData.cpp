//===============================================================================
//   Authors: Robert SCHMIEDER, Jeff SADURAL and Adrian CANTU
//     Computational Science Research Center @ SDSU, CA
//
//   File: prinseq-lite
//   Date: 2012-05-28
//   Version: 0.19.2 lite
//
//   Usage:
//      prinseq-lite [options]
//
//      Try 'prinseq-lite -h' for more information.
//
//   Purpose: PRINSEQ will help you to preprocess your genomic or metagenomic
//             sequence data in FASTA or FASTQ format. The lite version does not
//             require any non-core perl modules for processing.
//
//    Bugs: Please use http://sourceforge.net/tracker/?group_id=315449
//
//===============================================================================

#include "SequenceData.h"

SequenceData::SequenceData(){
    
}
SequenceData::SequenceData(string label,string sequence){
    sequenceID = label;
    this->sequence = sequence;
}

SequenceData::SequenceData(string label,string sequence, bool amino){
    
}

void SequenceData::SetSeqClass(string label,string sequence){
    
}
void SequenceData::SetID(string ID){
    this->sequenceID = ID;
}
void SequenceData::SetDNA(string sequence){
    this->sequence = sequence;
}
void SequenceData::SetSeqLength(int sequenceLength){
    this->sequenceLength = sequenceLength;
}
void SequenceData::SetAmino(bool amino){
    
}

string SequenceData::GetID(){
    return sequenceID;
}

string SequenceData::GetDNASeq(){
    return sequence;
}

long SequenceData::GetSeqLength(){
    return sequence.length();
}

bool SequenceData::GetAmino(){
    return amino;
}

void SequenceData::TrimSeqLeft(int trimValue){
    sequence.erase(sequence.begin(), sequence.begin() + trimValue);
}
void SequenceData::TrimSeqRight(int trimValue){
    sequence.erase(sequence.end(), sequence.end() - trimValue);
}
