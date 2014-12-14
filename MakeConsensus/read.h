#ifndef READ_H
#define READ_H
#include "qstring.h"

class Read
{
private:
    int length;
    int offset;
    QString sequence;
public:
    Read();
    void setLength(int n);
    void setOffset(int n);
    void setSequence(QString seq);
    int getLength();
    int getOffset();
    QString getSequence();


};

#endif // READ_H
