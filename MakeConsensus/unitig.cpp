#include "unitig.h"
#include <string>

using namespace std;


Unitig::Unitig()
{
}

void Unitig::setConsensus(string consensus){
    this->consensus=consensus;
}

string Unitig::getConsensus(){
    return consensus;
}
