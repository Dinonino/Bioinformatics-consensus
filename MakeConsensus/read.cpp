#include "read.h"

void Read::setLength(int n){
    length=n;
}

void Read::setOffset(int n){
    offset=n;
}

void Read::setSequence(QString seq){
    sequence=seq;
}

int Read::getLength(){
    return length;
}

int Read::getOffset(){
    return offset;
}

QString Read::getSequence(){
    return sequence;
}

Read::Read(){

}
