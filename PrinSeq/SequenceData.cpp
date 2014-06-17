//
//  SequenceData.cpp
//  PrinSeq
//
//  Created by Jeffrey Sadural on 6/17/14.
//  Copyright (c) 2014 Jeffrey Sadural. All rights reserved.
//

#include "SequenceData.h"

SequenceData::SequenceData(){
    
}
SequenceData::SequenceData(string label,string sequence){
    
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
    
}
void SequenceData::TrimSeqRight(int trimValue){
    
}