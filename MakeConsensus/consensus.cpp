#include "consensus.h"
#include "nucleic_codes.h"
#include "iostream"


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
    return sequence.size();
}

int Consensus::getOffset(){
    return offset;
}

void Consensus::append(column col)
{
    sequence.push_back(col);
}

column Consensus::getColumn(int index)
{
    return sequence.at(index);
}

string Consensus::getSequence(){
    string seq="";
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

    /*cerr << offset << endl;
    cerr << length << endl;
    cerr << this->getLength() << endl;
    cerr <<this->getOffset();*/
    int i=0;
    for(i=0;i<length;i++) {
        if((i+offset)<(this->offset) || (i+offset)>=(this->offset+this->getLength())) {
            column c;
            c.chatAt=nc.all;
            c.total=0;
            newCons.append(c);
        } else {
            newCons.append(sequence.at(offset+i-this->offset));
        }
    }
    return newCons;

}

