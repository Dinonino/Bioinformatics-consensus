#include "consensus.h"
#include "nucleic_codes.h"

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

column Consensus::getColumn(int index)
{
    return sequence.at(index);
}

QString Consensus::getSequence(){
    QString seq="";
    int i;
    for(i=0;i<this->getLength();i++) {
        seq+=sequence.at(i).chatAt;
    }
    return seq;
}

Consensus Consensus::getSubConsensus(int offset, int length)
{
    Nucleic_codes nc;
    Consensus newCons;
    newCons.setOffset(offset);

    int i=0;
    for(i=0;i<length;i++) {
        if((i+offset)<this->offset || (i+offset)>=(this->offset+this->getLength())) {
            column c;
            c.chatAt=nc.all;
            newCons.append(c);//TODO append end dash
            continue;
        } else {
            newCons.append(sequence.at(offset+i));
        }
    }
    return newCons;

}

