#include <QCoreApplication>
#include <iostream>
#include <fstream>
#include <QString>
#include <QList>
#include <QFile>
#include <QTextStream>
#include <QStringList>
#include <QDebug>

#include "read.h"
#include "unitig.h"
#include "realigner.h"
#include "nucleic_codes.h"
#include "consensus.h"

int main()
{
    Realigner realigner;

    QList<Unitig> unitigs;

    QFile readsFile("C:/Users/Pickles/Desktop/Bioinformatiks/Bioinformatics-consensus/lib/reads.2k.10x.fasta");
    QFile layoutFile ("C:/Users/Pickles/Desktop/Bioinformatiks/Bioinformatics-consensus/lib/layouts.afg");

    QList<QString> readsStringList;
    QString readString="";


    if (!readsFile.open(QIODevice::ReadOnly | QIODevice::Text)){
        std::cout << "Unable to open file";
        exit(-1);
    }

    QTextStream inReads(&readsFile);

    while (!inReads.atEnd()) {

        QString line = inReads.readLine();

        while(line.at(0)!='>'){

            readString+=line;
            line = inReads.readLine();

            if(inReads.atEnd()) break;
        }

        if(readString.isEmpty()) continue;

        readsStringList << readString;
        readString="";

    }

    readsFile.close();


    if (!layoutFile.open(QIODevice::ReadOnly | QIODevice::Text)){
        std::cout << "Unable to open file";
        exit(-1);
    }

    QTextStream inUnitigs(&layoutFile);



    while (!inUnitigs.atEnd()) {

        QString line = inUnitigs.readLine();
        if(line.compare("{LAY")==0) continue;
        QList<Read> unitigSequences;

            while(line.compare("}")!=0){
                Read* sequence=new Read();
                line=inUnitigs.readLine();
                QString clr=line.split(":")[1];
                QStringList startEnd=clr.split(",");
                int start=startEnd[0].toInt();
                int end=startEnd[1].toInt();

                line=inUnitigs.readLine();
                int offset=line.split(":")[1].toInt();

                line=inUnitigs.readLine();
                int src=line.split(":")[1].toInt();
                inUnitigs.readLine();

                line=inUnitigs.readLine();

                sequence->setOffset(offset);

                if(start<end){
                    sequence->setLength(end-start);
                    sequence->setSequence(readsStringList.at(src-1).mid(0,end-start));

                } else {
                    sequence->setLength(start-end);
                    std::string cstr = readsStringList.at(src-1).toStdString();
                    std::reverse(cstr.begin(), cstr.end());
                    QString reverseSequence=QString::fromStdString(cstr);
                    sequence->setSequence(reverseSequence.mid(0,start-end));

                }
                unitigSequences << *sequence;
                delete sequence;

            }

            Unitig unitig;
            unitig.sequences = unitigSequences;
            unitig.setStartEnd();
            unitigs << unitig;


    }

    layoutFile.close();


    int i=10; //npr. i=10
    int score;
    for(Unitig unitig : unitigs){
        QString consensusA=realigner.getAndScoreConsensus(unitig,&score);

        //ova funkcija je samo za ispis:
        for(int k=0; k<i && k<unitig.sequences.size(); k++){
            Read sequence=unitig.sequences.at(k);
            QString seq=sequence.getSequence();
            for(int z=0; z<sequence.getOffset();z++){
                seq.prepend(' ');
            }
            qDebug() << k+1 << ". "<<seq;

        }

         qDebug() << consensusA  << " Consensus A ";

        for(int k=0; k<i && k<unitig.sequences.size(); k++){
            Read sequence=unitig.sequences.at(k);
            unitig.removeSequence(k);

            Consensus consensusB=realigner.getConsensus2(unitig);
            consensusB.setOffset(unitig.getStart());

            sequence=realigner.align(consensusB, sequence, 0.1);

            qDebug() << "sequence " << k+1 << ".";

            unitig.insertSequnce(k, sequence);

            int newScore;
           qDebug() << consensusB.getSequence() <<" Consensus B";
            consensusA=realigner.getConsensus2(unitig).getSequence();

            //for petlja samo za ispis:
             qDebug() << "New alignment:";
            for(int j=0; j<i && j<unitig.sequences.size(); j++){
                Read sequence=unitig.sequences.at(j);
                QString seq=sequence.getSequence();
                for(int i=0; i<sequence.getOffset();i++){
                    seq.prepend(' ');
                }
                if(sequence.getOffset()>=0 && unitig.getStart()<0){
                    for(int i=unitig.getStart();i<0;i++){
                        seq.prepend('.');
                    }
                } else if(sequence.getOffset()<0 && sequence.getOffset()!=unitig.getStart()){
                    for(int i=0; i<sequence.getOffset()-unitig.getStart(); i++){
                        seq.prepend('.');
                    }
                }
                qDebug() << j+1 << "."<<seq;

            }

            qDebug()  << consensusA<< " New consensus A, offset: " << unitig.getStart();
           // if (newScore > score) break;   todo:change scoring function

        }

    }


    return 0;
}

