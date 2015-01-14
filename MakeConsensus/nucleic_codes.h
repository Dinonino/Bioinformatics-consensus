#ifndef NUCLEIC_CODES_H
#define NUCLEIC_CODES_H
#include <map>

using namespace std;

class Nucleic_codes
{
public:
    Nucleic_codes();
    map<char,unsigned char> getCharToByteMap();
    map<unsigned char, char> getByteToCharMap();
    char getCharFromByte(unsigned char byte);
    unsigned char getByteFromChar(char character);

    static const unsigned char dash;
    static const unsigned char A;
    static const unsigned char C;
    static const unsigned char G;
    static const unsigned char T;
    static const unsigned char all;
private:
    map<char, unsigned char> charToByteMap;
    map<unsigned char, char> byteToCharMap;
    void fillCharToByteMap();
    void fillByteToCharMap();
};



#endif // NUCLEIC_CODES_H
