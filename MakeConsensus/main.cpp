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
    string readsLocation = "C:/Users/Josipa/Desktop/gitprojekti/Bioinformatics-consensus/lib/reads.2k.10x.fasta\0";
    string layoutLocation = "C:/Users/Josipa/Desktop/gitprojekti/Bioinformatics-consensus/lib/layouts.afg\0";



    vector<string> readsStringList;
    string readString="";

    readsFile.open(readsLocation.c_str());

    if (!readsFile.is_open()){
        std::cout << "Unable to open file";
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
   // std::cout << readsStringList.size();

    layoutFile.open(layoutLocation.c_str());
    if (!layoutFile.is_open()){
        std::cout << "Unable to open file";
        exit(-1);
    }





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
            unitigs.push_back(unitig);

    }

    layoutFile.close();


    int i=10;
    string consensusA;
    Unitig unitig;
    for(int i=0; i< unitigs.size(); i++){
    	unitig = unitigs.at(i);
        consensusA=realigner.getConsensus(unitig);
        cout << "Consensus A : " << consensusA;
        // TODO: score consenusA
        for(int k=0; k<i && k<unitig.sequences.size(); k++){ // or until score consensusA increase
            Read sequence=unitig.sequences.at(k);

            unitig.sequences.erase(unitig.sequences.begin() + k);
            string consensusB=realigner.getConsensus(unitig);
            cout << "Consensus B : " << consensusB;
           // TODO : align k sequence with consensuB : sequence=realigner.align(sequence, consensusB);
            unitig.sequences.insert(unitig.sequences.begin() + k, sequence);
            consensusA=realigner.getConsensus(unitig);
            //TODO: score consensusA
        }
    }


    return 0;
}
