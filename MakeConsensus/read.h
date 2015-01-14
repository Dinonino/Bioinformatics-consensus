#ifndef READ_H
#define READ_H

#include <string>

using namespace std;

/**
 * @brief The Read class Klasa koja predstavlja jedno očitanje unutar unitiga.
 */
class Read
{
private:
    int length;
    int offset;
    string sequence;
public:
    Read();
    void setLength(int n);
    void setOffset(int n);
    void setSequence(string seq);
    int getLength();
    int getOffset();
    string getSequence();


};

#endif // READ_H
