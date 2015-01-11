#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <list>
#include <vector>

#include "read.h"
#include "unitig.h"
#include "realigner.h"
#include "nucleic_codes.h"
#include "consensus.h"

using namespace std;

 vector<std::string> &split(const string &s, char delim, vector<string> &elems) {
        stringstream ss(s);
        string item;
        while (getline(ss, item, delim)) {
            if(!item.empty())
                elems.push_back(item);
        }
        return elems;
    }

    vector<std::string> split(const string &s, char delim) {
        vector<string> elems;
        split(s, delim, elems);
        return elems;
    }

int main()
{
    Realigner realigner;

    vector<Unitig> unitigs;

    ifstream readsFile;
    ifstream layoutFile;
    string readsLocation = "/home/dino/Documents/Bioinformatics-consensus/lib/reads.2k.10x2.fasta";
    string layoutLocation = "/home/dino/Documents/Bioinformatics-consensus/lib/layouts2.afg";

    //string readsLocation = "C:/Users/Josipa/Desktop/gitprojekti/Bioinformatics-consensusMasterWorking/lib/reads.2k.10x2.fasta";
    //string layoutLocation = "C:/Users/Josipa/Desktop/gitprojekti/Bioinformatics-consensusMasterWorking/lib/layouts2.afg";

    vector<string> readsStringList;
    string readString="";

    readsFile.open(readsLocation.c_str(), ifstream::in);
    layoutFile.open(layoutLocation.c_str(), ifstream::in);



    if (!readsFile.is_open()){
        std::cout << "Unable to open reaads file";
        exit(-1);
    }

    if (!layoutFile.is_open()){
        std::cout << "Unable to open layout file";
        exit(-1);
    }


    string line;
    while (!readsFile.eof()) {
        getline(readsFile, line);
        while(line.at(0)!='>'){
            readString+=line;
            getline(readsFile, line);
            if(readsFile.eof()) break;
        }
        if(readString.length() == 0) continue;
        readsStringList.push_back(readString);
        readString="";
    }

    readsFile.close();


    while (!layoutFile.eof()) {
        getline(layoutFile, line);
        if(line.compare("{LAY")==0) continue;
        vector<Read> unitigSequences;
        while(line.compare("}")!=0){
            Read* sequence=new Read();
            getline(layoutFile, line);
            string clr=split(line, ':')[1];
            vector<string> startEnd=split(clr, ',');
            int start=atoi(startEnd[0].c_str());
            int end=atoi(startEnd[1].c_str());
            getline(layoutFile, line);
            int offset= atoi(split(line, ':')[1].c_str());
            getline(layoutFile, line);
            int src=atoi(split(line, ':')[1].c_str());;
            getline(layoutFile, line);
            getline(layoutFile, line);
            sequence->setOffset(offset);
            if(start<end){
                sequence->setLength(end-start);
                sequence->setSequence(readsStringList[src-1].substr(0,end-start));
            } else {
                sequence->setLength(start-end);
                string cstr = readsStringList[src-1];
                reverse(cstr.begin(), cstr.end());
                string reverseSequence=cstr;
                sequence->setSequence(reverseSequence.substr(0,start-end));
            }
            // qDebug() << sequence->getSequence();
            unitigSequences.push_back(*sequence);
            delete sequence;
            }
        Unitig unitig;
        unitig.sequences = unitigSequences;

        for(int t=0; t<unitig.sequences.size(); t++){
            cout << unitig.sequences.at(t).getSequence() << endl;
        }

        unitig.setStartEnd();
        unitigs.push_back(unitig);
    }

    layoutFile.close();


    int iterationNumber=10; //npr. i=10, još ga ne koristimo dolje zasad
    int score;
    Unitig unitig;
    string consensusA;
    for(int i=0; i < unitigs.size() ; i++){

        unitig = unitigs.at(i);

        consensusA = realigner.getAndScoreConsensus(unitig,&score);
        Read sequence;
        string seq;
        //ova funkcija je samo za ispis:
        for(int k=0; k<unitig.sequences.size(); k++){
            sequence=unitig.sequences.at(k);
            seq=sequence.getSequence();
            for(int z=0; z<sequence.getOffset();z++){
                seq.insert(0, " ");
            }
           cout << k+1 << ". " << seq << endl;

        }

        std::cout << consensusA  << " Consensus A ";

        for(int k=0; k<=i && k<unitig.sequences.size(); k++){
            sequence=unitig.sequences.at(k);
            unitig.removeSequence(k);

            Consensus consensusB=realigner.getConsensus2(unitig);
            consensusB.setOffset(unitig.getStart());

            sequence=realigner.align(consensusB, sequence, 0.1);

            std::cout << "sequence " << k+1 << "." << endl;

            unitig.insertSequnce(k, sequence);

            int newScore;
            std::cout << consensusB.getSequence() <<" Consensus B" << endl;
            consensusA=realigner.getConsensus2(unitig).getSequence();

            //for petlja samo za ispis:
            std::cout << "New alignment:" << endl;
            for(int j=0; j<i && j<unitig.sequences.size(); j++){
                sequence=unitig.sequences.at(j);
                string seq=sequence.getSequence();
                for(int i=0; i<sequence.getOffset();i++){
                    seq.insert(0, "  ");
                }
                if(sequence.getOffset()>=0 && unitig.getStart()<0){
                    for(int i=unitig.getStart();i<0;i++){
                        seq.insert(0, ".");
                    }
                } else if(sequence.getOffset()<0 && sequence.getOffset()!=unitig.getStart()){
                    for(int i=0; i<sequence.getOffset()-unitig.getStart(); i++){
                        seq.insert(0, ".");
                    }
                }
                std::cout << j+1 << "." << seq << endl;

            }
            std::cout << consensusA << " New consensus A, offset: " << unitig.getStart();
           // if (newScore > score) break;   todo:change scoring function

        }

    }


    return 0;
}

