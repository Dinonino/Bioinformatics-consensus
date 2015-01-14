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
    string getAndScoreConsensus(Unitig unitig, double* score);
    Read align(Consensus consensusB, Read sequence, double E);
    Consensus getConsensus2(Unitig unitig);
private:
    char getAndScoreColumnConsensus(Unitig unitig, int colNum, int *score1, double *score2);
    void align(Read sequence,string consensus,int e);
    column getColumnConsensus2(Unitig unitig, int colNum);

};

#endif // REALIGNER_H
