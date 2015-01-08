#ifndef CONSENSUS_H
#define CONSENSUS_H
#include<QString>
#include<QList>
#include<QMap>

typedef struct {
    char chatAt;
    QMap<QChar,int> freq;
    int total;
} column;

class Consensus
{
private:
    QString stringSeq;
    int offset;
    QList<column> sequence;
public:
    Consensus();
    ~Consensus();
    void setOffset(int n);
    int getLength();
    int getOffset();
    void append(column col);
    column getColumn(int index);
    QString getSequence();
    Consensus getSubConsensus(int offset,int length);
};

#endif // CONSENSUS_H
