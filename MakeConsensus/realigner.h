#ifndef REALIGNER_H
#define REALIGNER_H

#include "unitig.h"
#include "consensus.h"
#include "read.h"

using namespace std;

class Realigner
{
public:
    Realigner();
    string getConsensus(Unitig unitig);
    string getAndScoreConsensus(Unitig unitig, int* score);
    Read align(Consensus consensusB, Read sequence, double E);
    Consensus getConsensus2(Unitig unitig);
private:
    char getColumnConsensus(Unitig unitig, int colNum);
    char getAndScoreColumnConsensus(Unitig unitig, int colNum, int *score);
    void align(Read sequence,string consensus,int e);
    column getColumnConsensus2(Unitig unitig, int colNum);

};

#endif // REALIGNER_H
