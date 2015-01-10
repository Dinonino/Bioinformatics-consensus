#include "unitig.h"

#include <string>

Unitig::Unitig()
{
}

void Unitig::setConsensus(string consensus){
    this->consensus=consensus;
}

string Unitig::getConsensus(){
    return consensus;
}

int Unitig::getEnd(){
    return end;
}

int Unitig::getStart(){
    return start;
}

void Unitig::removeSequence(int k){
    sequences.erase(sequences.begin() + k);

    setStartEnd();
}

void Unitig::insertSequnce(int k, Read sequence){
    sequences.insert(sequences.begin() + k,sequence);

    setStartEnd();
}

void Unitig::setStartEnd(){
    int endOfUnitig=0;
    int startOfUnitig=((Read)sequences.at(0)).getOffset();

    for (int i=0; i< sequences.size(); i++){
        if(sequences.at(i).getLength()+sequences.at(i).getOffset() > endOfUnitig){
            endOfUnitig=sequences.at(i).getLength()+sequences.at(i).getOffset();
        }
        if(sequences.at(i).getOffset()<startOfUnitig){
           startOfUnitig=sequences.at(i).getOffset();
        }
    }

    start=startOfUnitig;
    end=endOfUnitig;
}
