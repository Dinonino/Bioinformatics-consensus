#include "realigner.h"
#include "nucleic_codes.h"
#include <QMap>
#include <QMapIterator>
#include <QChar>

#include <QDebug>

Realigner::Realigner()
{
}

char Realigner::getColumnConsensus(Unitig unitig, int colNum){
    Nucleic_codes nc;
    QMap<QChar,int> frequency;
    for(Read sequence : unitig.sequences){
       // qDebug() << sequence.getSequence() << " "<< sequence.getOffset();
        if(sequence.getOffset()<= colNum && sequence.getOffset()+sequence.getLength() > colNum){
          //  qDebug() << sequence.getSequence().at(colNum-sequence.getOffset());
            frequency[sequence.getSequence().at(colNum-sequence.getOffset())]++;
          //  qDebug() << "OVAJ";
        }
    }

    int max=0;
    foreach(int value, frequency){
        if(value>max){
            max=value;
        }
    }

    if(max==0) return nc.getCharFromByte(nc.A | nc.C | nc.G | nc.T | nc.dash);

    // qDebug() << "\nmax: " << max;
    char result=0;
    QMapIterator<QChar, int> i(frequency);
    while (i.hasNext()) {
        i.next();
        if(i.value()==max){
            result= result | nc.getByteFromChar(i.key().toLatin1());

        }
    }
   // qDebug() << "------ " << nc.getCharFromByte(result) ;
    return nc.getCharFromByte(result);

}

QString Realigner::getConsensus(Unitig unitig){
    Nucleic_codes nc;
    int endOfUnitig=0;
    int startOfUnitig=((Read)unitig.sequences.at(0)).getOffset();
    for(Read sequence : unitig.sequences){
        if(sequence.getLength()+sequence.getOffset() > endOfUnitig){
            endOfUnitig=sequence.getLength()+sequence.getOffset();
        }
        if(sequence.getOffset()<startOfUnitig){
           startOfUnitig=sequence.getOffset();
        }
    }

    QString consensus="";
    for(int i=0; i<startOfUnitig ; i++){
        consensus.append(" ");  // or '-', i'm not sure (here consensus doesn't exist but we keep track of offset of consensus)
    }
    for(int i=startOfUnitig; i<endOfUnitig; i++){
        QString columnConsensus=getColumnConsensus(unitig, i);
        consensus.append(columnConsensus);
    }

    return consensus;
}
