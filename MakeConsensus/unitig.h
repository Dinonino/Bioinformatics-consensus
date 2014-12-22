#ifndef UNITIG_H
#define UNITIG_H
#include "qlist.h"
#include "read.h"
#include "qstring.h"

class Unitig
{
public:
    Unitig();
    QList<Read> sequences;
    QString getConsensus();
    void setConsensus(QString consensus);
    int getStart();
    int getEnd();
    void removeSequence(int k);
    void insertSequnce(int k, Read sequence);
    void setStartEnd();
private:
    QString consensus;
    int start;
    int end;

};

#endif // UNITIG_H
