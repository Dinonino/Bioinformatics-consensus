#ifndef UNITIG_H
#define UNITIG_H
#include "read.h"
#include <vector>


class Unitig
{
public:
    Unitig();
    vector<Read> sequences;
    string getConsensus();
    void setConsensus(string consensus);
    int getStart();
    int getEnd();
    void removeSequence(int k);
    void insertSequnce(int k, Read sequence);
    void setStartEnd();
private:
    string consensus;
    int start;
    int end;

};

#endif // UNITIG_H
