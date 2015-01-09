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
        ret.total=0;
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
    int diffLen=(int)(E*sequence.getLength());
    int len=sequence.getLength()+diffLen*2;
    double* nw=(double*)malloc(sizeof(double)*(sequence.getLength()+1)*(len+1));
    char* direction=(char*)malloc((sequence.getLength()+1)*(len+1));




    Consensus cons=consensusB.getSubConsensus(sequence.getOffset()-diffLen,len);
    cons.setOffset(consensusB.getOffset()-diffLen);
    int i,j;
    double match,insert,max=-1;
    int consLen=cons.getLength()+1;
    int seqLen=sequence.getLength()+1;
    for(i=0;i<consLen;i++) {
        nw[i]=0;
        direction[i]=-1;
    }
    for(i=0;i<seqLen;i++) {
        nw[i*seqLen]=1000000;
        direction[i*seqLen]=-1;
    }
    for(i=1;i<seqLen;i++) {
        for(j=1;j<consLen;j++) {
            column col=cons.getColumn(j-1);
            if(col.total==0) {
                double match,insert;
                match=nw[(i-1)*consLen+j-1];
                insert=nw[i*consLen+j-1]+0.25;
                if(match>insert) nw[i*consLen+j]=insert;
                else nw[i*consLen+j]=match;
            } else {
                char seqChar=sequence.getSequence().at(i-1).toLatin1();
                char  seqByte=nc.getByteFromChar(seqChar);
                char consByte=nc.getByteFromChar(col.chatAt);
                double matched=1,inserted=1;
                if((consByte&seqByte)>0) matched=0;
                if((consByte&dash)>0) inserted=0;
                match=nw[(i-1)*consLen+j-1]+0.5*matched+0.5*(1-col.freq[seqChar]/col.total);
                insert=nw[i*consLen+j-1]+0.25+0.5*inserted+0.5*(1-col.freq['-']/col.total);
                if(match>insert) {
                    nw[i*consLen+j]=insert;
                    direction[i*consLen+j]=0;
                } else {
                    nw[i*consLen+j]=match;
                    direction[i*consLen+j]=1;
                }
            }
        }
    }

    int maxCol;
    if(sequence.getOffset()-diffLen<consensusB.getOffset()) {//sufiks-prefiks

        int sufixLen=0;
        while(sufixLen!=cons.getLength()-1 || cons.getColumn(sufixLen).total==0) {
            sufixLen++;
        }
        int out=0;
        for(i=1;i<consLen;i++) {
            if(i-1<(sufixLen)) out=sufixLen-(i-1);
            if(max==-1 || nw[(seqLen-1)*consLen+i]<max) {
                max=nw[(seqLen-1)*consLen+i]/(seqLen-1-out);
                maxCol=i;
            }
            out=0;
        }

    } else if(sequence.getOffset()+len-diffLen>=consensusB.getOffset()+consensusB.getLength()) {//prefiks-sufiks
        int prefixLen=0;
        while(prefixLen!=cons.getLength()-1 || cons.getColumn(consLen-2-prefixLen).total==0) {
            prefixLen++;

        }
        int end=cons.getLength()-1-prefixLen;
        int out=0;
        for(i=1;i<consLen;i++) {
            if(end<(i-1)) out=i-1-end;
            if(max==-1 || nw[(seqLen-1)*consLen+i]<max) {
                max=nw[(seqLen-1)*consLen+i]/(seqLen-1-out);
                maxCol=i;
            }
            out=0;
        }

    } else {//podniz
       for(i=1;i<consLen;i++) {
           if(max==-1 || nw[(seqLen-1)*consLen+i]<max) {
               max=nw[(seqLen-1)*consLen+i];
               maxCol=i;
           }
       }
    }
    QString final="";
    i=seqLen-1;
    j=maxCol;
    QString seqString=sequence.getSequence();
    while(i>0) {
        if(direction[i*consLen+j]) {
           final.insert(0,seqString.at(i-1));
           i--;
        } else {
           final.insert(0,"-");
        }
        j--;
    }
    Read returnRead;
    returnRead.setSequence(final);
    returnRead.setOffset(consensusB.getOffset()-diffLen+j);
    free(direction);
    free(nw);

















}
