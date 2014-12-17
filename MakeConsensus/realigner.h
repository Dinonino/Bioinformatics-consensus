#ifndef REALIGNER_H
#define REALIGNER_H

#include "unitig.h"
#include "qstring.h"


class Realigner
{
public:
    Realigner();
    char getColumnConsensus(Unitig unitig, int colNum);
    QString getConsensus(Unitig unitig);
private:

};

#endif // REALIGNER_H
