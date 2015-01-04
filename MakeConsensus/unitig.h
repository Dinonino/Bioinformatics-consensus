#ifndef UNITIG_H
#define UNITIG_H

#include "read.h"
#include <vector>
#include <string>

using namespace std;

class Unitig
{
public:
    Unitig();
    vector<Read> sequences;
    string getConsensus();
    void setConsensus(string consensus);
private:
    string consensus;
};

#endif // UNITIG_H
