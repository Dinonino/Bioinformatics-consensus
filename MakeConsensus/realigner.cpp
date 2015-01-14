#include "realigner.h"
#include "nucleic_codes.h"

#include <map>
#include <stdlib.h>
#include <iostream>

using namespace std;

Realigner::Realigner()
{
}


char Realigner::getAndScoreColumnConsensus(Unitig unitig, int colNum, int *score1, double *score2){
    Nucleic_codes nc;
    map<char,int> frequency;
    Read sequence;
    int i;
    int numberOfColumns=0;
    for(i=0; i < unitig.sequences.size() ; i++){
        sequence = unitig.sequences.at(i);
        if(sequence.getOffset()<= colNum && sequence.getOffset()+sequence.getLength() > colNum){
            frequency[sequence.getSequence().at(colNum-sequence.getOffset())]++;
	    numberOfColumns++;	
        }
    }

    int max=0;
    map<char, int>::iterator iter;
    for (iter = frequency.begin(); iter != frequency.end(); iter++){
        if(iter->second>max){
            max=iter->second;
        }
    }
    *score1=0;
    *score2=0;	
    if(max==0) return nc.getCharFromByte(nc.A | nc.C | nc.G | nc.T | nc.dash);


    char result=0;

    for (iter = frequency.begin(); iter != frequency.end(); iter++){
        if(iter->second==max){
            result = result | nc.getByteFromChar(iter->first);

        } else {
            *score1=*score1+iter->second;
        }
    }
    *score2=*score1/numberOfColumns;	
    return nc.getCharFromByte(result);

}

Consensus Realigner::getConsensus2(Unitig unitig)
{
    Consensus consensus;
 /*   for(int i=0; i<unitig.getStart() ; i++){
        consensus.append(" ");  // or '-' or '&', i'm not sure (here consensus doesn't exist but we keep track of offset of consensus)
    }*/

    for(int i=unitig.getStart(); i<unitig.getEnd(); i++){
        column columnConsensus=getColumnConsensus2(unitig, i);
        consensus.append(columnConsensus);
    }

    return consensus;
}

column Realigner::getColumnConsensus2(Unitig unitig, int colNum)
{
    Nucleic_codes nc;
    column ret;
    map<char,int> frequency;
    frequency['g']=0;
    frequency['c']=0;
    frequency['a']=0;
    frequency['t']=0;
    frequency['-']=0;
    Read sequence;
    for(int i=0; i < unitig.sequences.size(); i++){
        sequence = unitig.sequences.at(i);
        if(sequence.getOffset()<= colNum && sequence.getOffset()+sequence.getLength() > colNum){
            frequency[sequence.getSequence().at(colNum-sequence.getOffset())]++;
        }
    }

    int max=0;
    int total=0;
    map<char, int>::iterator iter;
    for (iter = frequency.begin(); iter != frequency.end(); iter++){
        total+= iter->second;
        if(iter->second>max){
            max=iter->second;
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
    for (iter = frequency.begin(); iter != frequency.end(); iter++){
        if(iter->second==max){
            result= result | nc.getByteFromChar(iter->first);

        }
    }
   // qDebug() << "------ " << nc.getCharFromByte(result) ;
    ret.chatAt=nc.getCharFromByte(result);
    return ret;
}


string Realigner::getAndScoreConsensus(Unitig unitig, double *score){
    string consensus="";
 /*   for(int i=0; i<unitig.getStart() ; i++){
        consensus.append(" ");  // or '-' or '&', i'm not sure (here consensus doesn't exist but we keep track of offset of consensus)
    }*/

    int columnScore1;
    double columnScore2;
    *score=0;
    string columnConsensus;
    
    int score1=0;
    double score2=0;

    int i;
    for(i=unitig.getStart(); i<unitig.getEnd(); i++){

      //  ss << getAndScoreColumnConsensus(unitig, i,&columnScore);
      //  ss >> columnConsensus;
        consensus.append(1, getAndScoreColumnConsensus(unitig, i,&columnScore1, &columnScore2));
        score1=score1+columnScore1;
        score2=score2+columnScore2;
    }
    *score=0.5*score1+ 0.5*score2;

    return consensus;
}


/**
 * @brief Realigner::align Funkcija algin. Podešava trenutnu sekvencu s obzirom na ostale sekvence i njihov konsenzus. Trenutna sekvenca može poprimiti
 *  novi offset koji može biti udaljen za maksimalno E*seqence.length.
 * @param consensusB Konsenzus ostatka
 * @param sequence Sekvenca koja se ravna s ostatkom.
 * @param E Faktor maksimalnog pomaka sekvence, za najbolje rezultate predlažu se vrijednosti [0,0.1].
 * @return Poravnata sekvenca s ostalim sekvencama i njihovim konsenzusom.
 */
Read Realigner::align(Consensus consensusB, Read sequence, double E)
{
   // qDebug()<<"align funkcija";
    Nucleic_codes nc;
    char dash=nc.dash;
    int diffLen=(int)(E*sequence.getLength());
    int len=sequence.getLength()+diffLen*2;
    double* nw=(double*)malloc(sizeof(double)*(sequence.getLength()+1)*(len+1));
    char* direction=(char*)malloc((sequence.getLength()+1)*(len+1));


    Consensus cons=consensusB.getSubConsensus(sequence.getOffset()-diffLen,len);

    int i,j;
    double match,insert,max=-1;
    int consLen=cons.getLength()+1;
    int seqLen=sequence.getLength()+1;
    for(i=0;i<consLen;i++) {
        nw[i]=0;
        direction[i]=-1;
    }
    for(i=1;i<seqLen;i++) {
        nw[i*consLen]=100;
        direction[i*seqLen]=-1;
    }


    /* Poedešeni Needleman Wunsch*/
    /* Nije moguće umetati '-' u konsenzus,samo u sekvencu*/
    /* Prvi red je postavljen na 0, pošto ne tražimo najbolji položaj sekvence s obzirom na konsenzus*/
    /* Funkcija dobrote poprima vrijednosti [0,1] te dodatno kažnjava umetanje '-' s 0.25. Manji iznos je bolji niz*/
    /* U slučaju da sekvenca nije podniz konsenzusa(npr. preklapanje prefiksa ili sufiksa) match se nagrađuje s -0.25, a za najbolji algin se uzima
        onaj koji ima najmanji iznos funkcije dobrote po preklapajućem charu.*/
    /* Dodatno inicijaliziramo polje smjera popunjavanja matrice kako bismo se lakše vratili na početak. Naime, pošto su vrijednosti u matrici
        tipa double potrebno je zaokruživati vrijednosti da znamo otkrijemo put prema natrag. Zaokruživanje troši vrijeme, pa smo se odlučili na novu matricu
        smjera koja samo privremeno troši memoriju*/
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
                char seqChar=sequence.getSequence().at(i-1);
                char  seqByte=nc.getByteFromChar(seqChar);
                char consByte=nc.getByteFromChar(col.chatAt);
                double matched=1,inserted=1;
                //qDebug()<<seqChar<<","<<col.chatAt;
                if((consByte&seqByte)>0) matched=0;
                if((consByte&dash)>0) inserted=0;
                //qDebug()<<QString::number(matched);
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

    string str="";
    string number;
    /* ISPIS needleman wunsch matrice*/
    /*
    for(i=0;i<cons.getSequence().length();i++) {
        cout<<cons.getSequence().at(i)+"   ";
    }
    cout<<"\n";
    for(i=0;i<seqLen;i++) {
        for(j=0;j<consLen;j++) {
            convert << (nw[i*consLen+j]);
            number = convert.str();
            convert.str("");
            str+=number+" ";
            int k;
            for(k=number.length();k<6;k++) {
                str+=" ";
            }
        }
      cout << str<<"\n";
        str="";
    }*/




    int maxCol;

    /* Slučaj kada je moguće da je konsenzus ostatka podniz same sekvence, tj. sekvenca može biti dulja od konsenzusa.*/
    if(sequence.getOffset()-diffLen<consensusB.getOffset() && sequence.getOffset()+sequence.getLength()+diffLen>=consensusB.getOffset()+consensusB.getLength()) {
        int sufixLen=0,prefixLen=0;
        while(sufixLen<cons.getLength()-1 && cons.getColumn(sufixLen).total==0) {
            sufixLen++;
        }
        while(prefixLen<cons.getLength()-1 && cons.getColumn(consLen-2-prefixLen).total==0) {
            prefixLen++;

        }
        int end=cons.getLength()-1-prefixLen;
        int out=1;
        for(i=1;i<consLen;i++) {
            if((i-1)>sufixLen) out=i-sufixLen;
            if(end<(i-1)) out=seqLen-1-out;
            if(end<(i-1) && (i-seqLen)<sufixLen) out=consensusB.getLength();
            double current=nw[(seqLen-1)*consLen+i]/out-0.25*(out/(seqLen-1));
            if(max==-1 || current<max) {
                max=nw[(seqLen-1)*consLen+i]/out;
                maxCol=i;
            }
            out=1;
        }

    }
    /* Preklapanje sufiksa sekvence s prefiksom konsenzusa*/
    else if(sequence.getOffset()-diffLen<consensusB.getOffset()) {

        int sufixLen=0;
        while(sufixLen<cons.getLength()-1 && cons.getColumn(sufixLen).total==0) {
            sufixLen++;
        }
        int out=1;
        for(i=1;i<consLen;i++) {
            if((i-1)>sufixLen) out=i-sufixLen;
            if(out>=(seqLen-1)) out=seqLen-1;
            double current=nw[(seqLen-1)*consLen+i]/out-0.25*(out/(seqLen-1));
            if(max==-1 || current<max) {
                max=current;
                maxCol=i;
            }
            out=1;
        }

    } /* Preklapanje prefiksa sekvence s sufiksom konsenzusa*/
    else if(sequence.getOffset()+sequence.getLength()+diffLen>=consensusB.getOffset()+consensusB.getLength()) {
        int prefixLen=0;
        while(prefixLen<cons.getLength()-1 && cons.getColumn(consLen-2-prefixLen).total==0) {
            prefixLen++;

        }
        int end=cons.getLength()-1-prefixLen;
        int out=0;
        for(i=1;i<consLen;i++) {
            if(end<(i-1)) out=i-1-end;
            double current=nw[(seqLen-1)*consLen+i]/(seqLen-1-out)-0.25*(out/(seqLen-1));
            if(max==-1 || current<max) {
                max=current;
                maxCol=i;
            }
            out=0;
        }

    }
    /* Sekvenca je podniz lonsenzusa*/
    else {
       for(i=1;i<consLen;i++) {
           if(max==-1 || nw[(seqLen-1)*consLen+i]<max) {
               max=nw[(seqLen-1)*consLen+i];
               maxCol=i;
           }
       }
    }

    //trace maksimuma do pocetka
    string final="";
    i=seqLen-1;
    j=maxCol;


    string seqString=sequence.getSequence();
    while(i>0) {
        if(direction[i*consLen+j]) {
           str= seqString.at(i-1);
           final.insert(0, str);
           i--;
        } else {
           final.insert(0,"-");
        }
        j--;
    }
    int newOffset=cons.getOffset()+j;
    /*cerr << "novo" << endl;
    cerr << "maxCol: " << maxCol << endl;
    cerr << "j: " << j << endl;
    cerr << "newoffset: "<< newOffset << endl;
    cerr << "cons.offset: " << cons.getOffset() << endl;
    cerr << "difflen: "<< diffLen << endl;*/

    Read returnRead;
    returnRead.setSequence(final);
    returnRead.setOffset(newOffset);
    free(direction);
    free(nw);
    return returnRead;

















}
