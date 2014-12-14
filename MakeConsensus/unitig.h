#ifndef UNITIG_H
#define UNITIG_H
#include "qlist.h"
#include "read.h"

class Unitig
{
public:
    Unitig();
    QList<Read> sequences;
};

#endif // UNITIG_H
