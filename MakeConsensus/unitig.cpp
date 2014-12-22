#include "unitig.h"
#include <QString>

Unitig::Unitig()
{
}

void Unitig::setConsensus(QString consensus){
    this->consensus=consensus;
}

QString Unitig::getConsensus(){
    return consensus;
}

int Unitig::getEnd(){
    return end;
}

int Unitig::getStart(){
    return start;
}

void Unitig::removeSequence(int k){
    sequences.removeAt(k);

    setStartEnd();
}

void Unitig::insertSequnce(int k, Read sequence){
    sequences.insert(k,sequence);

    setStartEnd();
}

void Unitig::setStartEnd(){
    int endOfUnitig=0;
    int startOfUnitig=((Read)sequences.at(0)).getOffset();
    for(Read sequence : sequences){
        if(sequence.getLength()+sequence.getOffset() > endOfUnitig){
            endOfUnitig=sequence.getLength()+sequence.getOffset();
        }
        if(sequence.getOffset()<startOfUnitig){
           startOfUnitig=sequence.getOffset();
        }
    }

    start=startOfUnitig;
    end=endOfUnitig;
}
