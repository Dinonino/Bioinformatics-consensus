#ifndef CONSENSUS_H
#define CONSENSUS_H

#include <map>
#include <string>
#include <vector>

using namespace std;

typedef struct {
    char chatAt;
    map<char,int> freq;
    int total;
} column;

class Consensus
{
private:
    string stringSeq;
    int offset;
    vector<column> sequence;
public:
    Consensus();
    ~Consensus();
    void setOffset(int n);
    int getLength();
    int getOffset();
    void append(column col);
    column getColumn(int index);
    string getSequence();
    Consensus getSubConsensus(int offset,int length);
};

#endif // CONSENSUS_H
