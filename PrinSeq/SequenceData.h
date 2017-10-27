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

#ifndef __PrinSeq__SequenceData__
#define __PrinSeq__SequenceData__

#define VALID_BASES_AA "[ABCDEFGHIKLMNOPQRSTUVWYZXabcdefghiklmmopqrstuvwyzx*-]+"

#define VALID_BASES_NON_AA "[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+"

#include <iostream>
using namespace std;

class SequenceData{
public:
    
    SequenceData();
    SequenceData(string label,string sequence);
    SequenceData(string label,string sequence, bool amino);
    
    void SetSeqClass(string label,string sequence);
    void SetID(string sequenceID); // remove
    void SetDNA(string sequence); // remove
    void SetSeqLength(int sequenceLength); // remove
    void SetAmino(bool amino);
    
    string GetID();
    string GetDNASeq();
    long GetSeqLength();
    bool GetAmino();
    
    void TrimSeqLeft(int trimValue);
    void TrimSeqRight(int trimValue);
    
private:
    string sequenceID;
    string sequence;
    long sequenceLength;
    bool amino;
    
    
};
#endif /* defined(__PrinSeq__SequenceData__) */
