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
        //qDebug() << sequence.getSequence() << " "<< sequence.getOffset();
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

char Realigner::getAndScoreColumnConsensus(Unitig unitig, int colNum, int *score){
    Nucleic_codes nc;
    QMap<QChar,int> frequency;
    for(Read sequence : unitig.sequences){
        if(sequence.getOffset()<= colNum && sequence.getOffset()+sequence.getLength() > colNum){
            frequency[sequence.getSequence().at(colNum-sequence.getOffset())]++;
        }
    }

    int max=0;
    foreach(int value, frequency){
        if(value>max){
            max=value;
        }
    }
    *score=0;
    if(max==0) return nc.getCharFromByte(nc.A | nc.C | nc.G | nc.T | nc.dash);


    char result=0;

    QMapIterator<QChar, int> i(frequency);
    while (i.hasNext()) {
        i.next();
        if(i.value()==max){
            result= result | nc.getByteFromChar(i.key().toLatin1());

        } else {
            *score=*score+i.value();
        }
    }

    return nc.getCharFromByte(result);

}

QString Realigner::getConsensus(Unitig unitig){

    QString consensus="";
 /*   for(int i=0; i<unitig.getStart() ; i++){
        consensus.append(" ");  // or '-' or '&', i'm not sure (here consensus doesn't exist but we keep track of offset of consensus)
    }*/
    for(int i=unitig.getStart(); i<unitig.getEnd(); i++){
        QString columnConsensus=getColumnConsensus(unitig, i);
        consensus.append(columnConsensus);
    }

    return consensus;
}

QString Realigner::getAndScoreConsensus(Unitig unitig, int *score){
    QString consensus="";
 /*   for(int i=0; i<unitig.getStart() ; i++){
        consensus.append(" ");  // or '-' or '&', i'm not sure (here consensus doesn't exist but we keep track of offset of consensus)
    }*/

    int columnScore;
    *score=0;
    for(int i=unitig.getStart(); i<unitig.getEnd(); i++){
        QString columnConsensus=getAndScoreColumnConsensus(unitig, i,&columnScore);
        consensus.append(columnConsensus);
        *score=*score+columnScore;
    }

    return consensus;
}
