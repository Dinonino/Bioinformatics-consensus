#include "read.h"
#include <string>

void Read::setLength(int n){
    length=n;
}

void Read::setOffset(int n){
    offset=n;
}

void Read::setSequence(string seq){
    sequence=seq;
}

int Read::getLength(){
    return length;
}

int Read::getOffset(){
    return offset;
}

string Read::getSequence(){
    return sequence;
}

Read::Read(){
	length = 0;
	offset = 0;
}
