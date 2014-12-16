#ifndef REALIGNER_H
#define REALIGNER_H

#include "unitig.h"


class Realigner
{
public:
    Realigner();
    char getColumnConsensus(Unitig unitig, int colNum);
private:

};

#endif // REALIGNER_H
