#ifndef REALIGNER_H
#define REALIGNER_H

#include "unitig.h"
#include <string>

using namespace std;


class Realigner
{
public:
    Realigner();
    char getColumnConsensus(Unitig unitig, int colNum);
    string getConsensus(Unitig unitig);
private:

};

#endif // REALIGNER_H
