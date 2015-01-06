#include "consensus.h"

Consensus::Consensus()
{

}

Consensus::~Consensus()
{

}


void Consensus::setOffset(int n){
    offset=n;
}


int Consensus::getLength(){
    return sequence.length();
}

int Consensus::getOffset(){
    return offset;
}

void Consensus::append(column col)
{
    sequence.append(col);
}

QString Consensus::getSequence(){
    return stringSeq;
}

