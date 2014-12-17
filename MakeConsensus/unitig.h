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
private:
    QString consensus;
};

#endif // UNITIG_H
