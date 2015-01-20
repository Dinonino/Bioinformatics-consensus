#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <vector>
#include "read.h"
#include "unitig.h"
#include "realigner.h"
#include "nucleic_codes.h"
#include "consensus.h"
#include <algorithm>

using namespace std;

/**
  * @brief split Method for spliting string by delimiter
  * @param s String that needs to be split
  * @param delim Delimiter
  * @param elems Result of spliting string
  * @return
  */
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

/**
 * @brief writeConsensus Method for writing consensus in file and on standard output.
 * @param consensus Consensus
 * @param file
 */
void writeConsensus(Consensus consensus, ofstream* file){
    *file << "<\n";
    Nucleic_codes nc;

    string consensusSeq=consensus.getSequence();
    for(int i=0; i<consensusSeq.size(); i++){
        vector<string> bases;
        char consBase=consensusSeq.at(i);

        if(consBase=='N'){
            *file <<"N";
            continue;
        }

        unsigned char byte=nc.getByteFromChar(consBase);

        if((byte & nc.A) > 0){
            bases.push_back("A");
        }
        if((byte & nc.C) > 0){
            bases.push_back("C");
        }
        if((byte & nc.T) > 0){
            bases.push_back("T");
        }
        if((byte & nc.G)> 0){
            bases.push_back("G");
        }
        if((byte & nc.dash) > 0){
            bases.push_back("-");
        }

        if(bases.size()>1){
            *file <<"[";

            for(int j=0; j<bases.size(); j++){
                *file<<bases[j];
            }

            *file <<"]";
        }else {
            *file<<consBase;
        }


    }
    *file <<"\n";
}

/**
 * @brief printAlignment Method for printing unitig on standard output
 * @param unitig Unitig that will be printed
 */
void printAlignment (Unitig unitig){
    Read sequence;
    int j;
    for(j=0; j<unitig.sequences.size(); j++){
        sequence=unitig.sequences.at(j);
        string seq=sequence.getSequence();
        for(int i=0; i<sequence.getOffset();i++){
            seq.insert(0, " ");
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
        cout << j+1 << ". " << seq << endl;

    }

    cout << "--------------------------\n";
}

/**
 * @brief readReads Method used for reading reads from file
 * @param readsLocation Reads location path
 * @return Vector of reads
 */
vector<string> readReads(string readsLocation){
     ifstream readsFile;
     string readString="";
     vector<string> readsStringList;

     readsFile.open(readsLocation.c_str(), ifstream::in);

     if (!readsFile.is_open()){
         cout << "Unable to open reads file (first argument)";
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
     return readsStringList;
}

/**
 * @brief readLayout Method used for reading layout from file
 * @param layoutLocation Layout location path
 * @param readsStringList Vector of ureads
 * @return Vector of unitigs
 */
vector<Unitig> readLayout(string layoutLocation, vector<string> readsStringList){
    ifstream layoutFile;
    layoutFile.open(layoutLocation.c_str(), ifstream::in);
    vector<Unitig> unitigs;

    if (!layoutFile.is_open()){
        cout << "Unable to open layout file (second argument)";
        exit(-1);
    }

    string line;
    while (!layoutFile.eof()) {
        getline(layoutFile, line);
        if(line.empty()) break;
        if(line.compare("{LAY")==0) continue;
        vector<Read> unitigSequences;
        while(line.compare("}")!=0 ){
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
            unitigSequences.push_back(*sequence);
            delete sequence;

            }
        Unitig unitig;
        unitig.sequences = unitigSequences;


        unitig.setStartEnd();
        unitigs.push_back(unitig);
    }

    layoutFile.close();
    return unitigs;
}



int main(int argc, char* argv[])
{
    Realigner realigner;
    vector<string> readsStringList;
    vector<Unitig> unitigs;

    if(argc < 5){
        cout << "Not enough parameters. You must enter file path to reads, file path to layout, number of iterations and error rate";
        exit(-1);
    }

    int iterationNumber;
    try {
        iterationNumber=stoi(argv[3]);
    } catch(exception e){
        cout << "Number of iterations ( third argument) must be int number";
        exit(-1);
    }

    double errorRate;
    try {
        errorRate=stod(argv[4]);
    } catch(exception e){
        cout << "Error rate ( fourth argument) must be some real number";
        exit(-1);
    }


    ofstream consensusFile;
    ofstream oldConsensusFile;

    string readsLocation = argv[1];
    string layoutLocation = argv[2];
    string consensusResultLocation="consensus.txt";
    string oldConsensusLocation="old_consesnsus.txt";

    readsStringList=readReads(readsLocation);
    unitigs=readLayout(layoutLocation, readsStringList);

    consensusFile.open(consensusResultLocation.c_str(), ofstream::out);
    oldConsensusFile.open(oldConsensusLocation.c_str(), ofstream::out);



    Unitig unitig;
    Consensus consensusA;
    Read sequence;

    /*
     * Pseudokod:
     *
     za svaki unitig radi
        za k=1 do  k=broj koraka
            izdvoji k-to očitanje iz unitiga
            od ostalih očitanja napravi novi (privremeni) konsenzus
            poravnaj izdvojeno očitanje s privremenim konsezusom
            vrati izdvojeno očitanje u unitig
        izračunaj novi konsenzus unitiga i spremi ga

     */


    for(int i=0; i < unitigs.size() ; i++){

        unitig = unitigs.at(i);
        consensusA = realigner.getConsensus2(unitig);

        printAlignment(unitig);
        cout << "   " << consensusA.getSequence()  << " Consensus A\n ";
        writeConsensus(realigner.getConsensus2(unitig), &oldConsensusFile);

        for(int k=0; k<unitig.sequences.size() && k<=iterationNumber; k++){
            sequence=unitig.sequences.at(k);
            unitig.removeSequence(k);

            Consensus consensusB=realigner.getConsensus2(unitig);
            consensusB.setOffset(unitig.getStart());

            sequence=realigner.align(consensusB, sequence, errorRate);

            unitig.insertSequnce(k, sequence);

            consensusA=realigner.getConsensus2(unitig);

            cout << "\nNew alignment:" << endl;
            printAlignment(unitig);
            cout << "   " << consensusA.getSequence() << " New consensus A\n" ;

        }

        writeConsensus(consensusA, &consensusFile);
        cout<<"\nend of unitig\n\n";


    }

    consensusFile.close();
    oldConsensusFile.close();
    return 0;
}

