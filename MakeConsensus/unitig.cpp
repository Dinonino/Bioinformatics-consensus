#include "unitig.h"
#include <QString>

Unitig::Unitig()
{
}

void Unitig::setConsensus(QString consensus){
    this->consensus=consensus;
}

QString Unitig::getConsensus(){
    return consensus;
}
