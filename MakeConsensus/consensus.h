#ifndef CONSENSUS_H
#define CONSENSUS_H
#include<QString>
#include<QList>

typedef struct {
    QString chatAt;
    int mismatches;
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
    QString getSequence();
};

#endif // CONSENSUS_H
