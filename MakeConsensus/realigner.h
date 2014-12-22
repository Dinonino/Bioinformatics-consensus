#ifndef REALIGNER_H
#define REALIGNER_H

#include "unitig.h"
#include "qstring.h"


class Realigner
{
public:
    Realigner();
    QString getConsensus(Unitig unitig);
    QString getAndScoreConsensus(Unitig unitig, int* score);
private:
    char getColumnConsensus(Unitig unitig, int colNum);
    char getAndScoreColumnConsensus(Unitig unitig, int colNum, int *score);

};

#endif // REALIGNER_H
