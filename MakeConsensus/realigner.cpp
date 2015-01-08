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
    frequency['g']=0;
    frequency['c']=0;
    frequency['a']=0;
    frequency['t']=0;
    frequency['-']=0;
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

Consensus Realigner::getConsensus2(Unitig unitig)
{
    Consensus *consensus=new Consensus();
 /*   for(int i=0; i<unitig.getStart() ; i++){
        consensus.append(" ");  // or '-' or '&', i'm not sure (here consensus doesn't exist but we keep track of offset of consensus)
    }*/
    for(int i=unitig.getStart(); i<unitig.getEnd(); i++){
        column columnConsensus=getColumnConsensus2(unitig, i);
        consensus->append(columnConsensus);
    }

    return *consensus;
}

column Realigner::getColumnConsensus2(Unitig unitig, int colNum)
{
    Nucleic_codes nc;
    column ret;
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
    int total=0;
    foreach(int value, frequency){
        total+=value;
        if(value>max){
            max=value;
        }
    }
    /*int mismatch;
    foreach(int value, frequency){
        if(value<max){
            mismatch+=value;
        }
    }*/
    ret.total=total;
    ret.freq=frequency;

    if(max==0) {
        ret.chatAt=nc.getCharFromByte(nc.A | nc.C | nc.G | nc.T | nc.dash);
        return ret;
    }

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
    ret.chatAt=nc.getCharFromByte(result);
    return ret;
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

void Realigner::align(Consensus consensusB, Read sequence, double E)
{
    Nucleic_codes nc;
    char dash=nc.dash;
    char empty=nc.getCharFromByte(nc.A | nc.C | nc.G | nc.T | nc.dash);
    int len=sequence.getLength()+((int)(E*sequence.getLength()))*2;
    double* nw=(double*)malloc(sizeof(double)*(sequence.getLength()+1)*(len+1));
    if(sequence.getOffset()-sequence.getLength()*E<consensusB.getOffset()) {//sufiks-prefiks

    } else if(sequence.getOffset()+len-((len-sequence.getLength())/2)>=consensusB.getOffset()+consensusB.getLength()) {//prefiks-sufiks

    } else {//podniz
        Consensus cons=consensusB.getSubConsensus(sequence.getOffset()-sequence.getLength()*E,len);

        int i,j;
        double match,insert;
        int consLen=cons.getLength()+1;
        int seqLen=sequence.getLength()+1;
        for(i=0;i<consLen;i++) nw[i]=0;
        for(i=0;i<seqLen;i++) nw[i*seqLen]=1000000;
        for(i=1;i<seqLen;i++) {
            for(j=1;j<consLen;j++) {
                column col=cons.getColumn(j-1);
                char seqChar=sequence.getSequence().at(i-1).toLatin1();
                char  seqByte=nc.getByteFromChar(seqChar);
                char consByte=col.chatAt;
                double matched=0,inserted=0;
                if((consByte&seqChar)>0) matched=1;
                if((consByte&dash)>0) inserted=1;
                match=nw[(i-1)*seqLen+j-1]+0.5*matched+0.5*(col.freq[seqChar]/col.total);
                insert=nw[i*seqLen+j-1]+0.25+0.5*inserted+0.5*(col.freq['-']/col.total);
                if(match<insert) {
                    nw[i*seqLen+j]=insert;
                } else {
                    nw[i*seqLen+j]=match;
                }
            }
        }
    }
    free(nw);

















}
