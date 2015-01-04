#include "realigner.h"
#include "nucleic_codes.h"
#include <map>
#include <iterator>
#include <string>

using namespace std;

Realigner::Realigner()
{
}

char Realigner::getColumnConsensus(Unitig unitig, int colNum){
    Nucleic_codes nc;
    map<char,int> frequency;
    Read sequence;
    for(int i=0; i<unitig.sequences.size(); i++){
    	sequence = unitig.sequences[i];
       // qDebug() << sequence.getSequence() << " "<< sequence.getOffset();
        if(sequence.getOffset()<= colNum && sequence.getOffset()+sequence.getLength() > colNum){
          //  qDebug() << sequence.getSequence().at(colNum-sequence.getOffset());
            frequency[sequence.getSequence()[colNum-sequence.getOffset()]]++;
          //  qDebug() << "OVAJ";
        }
    }

    int max=0;
    map<char, int>::iterator iter;

    for (iter = frequency.begin(); iter != frequency.end(); ++iter){
    	if(iter->second>max){
    		max=iter->second;
    	}
    }



    if(max==0) return nc.getCharFromByte(nc.A | nc.C | nc.G | nc.T | nc.dash);


    char result=0;


    for (iter = frequency.begin(); iter != frequency.end(); ++iter){
    	if(iter->second==max){
    		result = result | nc.getByteFromChar(iter->first);
    	}
    }
    /*
    while (i.hasNext()) {
        i.next();
        if(i.value()==max){
            result= result | nc.getByteFromChar(i.key().toLatin1());

        }
    }
    */
   // qDebug() << "------ " << nc.getCharFromByte(result) ;
    return nc.getCharFromByte(result);

}

string Realigner::getConsensus(Unitig unitig){
    Nucleic_codes nc;
    int endOfUnitig=0;
    int startOfUnitig=((Read)unitig.sequences.at(0)).getOffset();
    Read sequence;
    for(int i=0; i<unitig.sequences.size(); i++){
    	sequence = unitig.sequences[i];
        if(sequence.getLength()+sequence.getOffset() > endOfUnitig){
            endOfUnitig=sequence.getLength()+sequence.getOffset();
        }
        if(sequence.getOffset()<startOfUnitig){
           startOfUnitig=sequence.getOffset();
        }
    }

    string consensus="";
    for(int i=0; i<startOfUnitig ; i++){
        consensus.append(" ");  // or '-', i'm not sure (here consensus doesn't exist but we keep track of offset of consensus)
    }
    for(int i=startOfUnitig; i<endOfUnitig; i++){
        consensus.push_back(getColumnConsensus(unitig, i));
    }

    return consensus;
}
