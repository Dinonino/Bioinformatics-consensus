#include "realigner.h"
#include "nucleic_codes.h"
#include <QMap>
#include <QMapIterator>
#include <QChar>

#include <QDebug>
#include <iostream>

Realigner::Realigner()
{


}

char Realigner::getColumnConsensus(Unitig unitig, int colNum){
    Nucleic_codes nc;

    QMap<QChar,int> frequency;
    for(Read sequence : unitig.sequences){
        //qDebug() << sequence.getSequence() << " "<< sequence.getOffset();
        if(sequence.getOffset()<= colNum && sequence.getOffset()+sequence.getLength() > colNum){
            frequency[sequence.getSequence().at(colNum-sequence.getOffset())]++;
            //qDebug() << "OVAJ";
        }
    }

    int max=0;
    foreach(int value, frequency){
        if(value>max){
            max=value;
        }
    }

     //qDebug() << "\nmax: " << max;
    char result=0;
    QMapIterator<QChar, int> i(frequency);
    while (i.hasNext()) {
        i.next();
        if(i.value()==max){
            result= result | nc.getByteFromChar(i.key().toLatin1());
            qDebug() << i.key();
        }
    }

    return nc.getCharFromByte(result);

}
