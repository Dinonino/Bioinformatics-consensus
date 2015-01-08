#ifndef NUCLEIC_CODES_H
#define NUCLEIC_CODES_H

#include  "qmap.h"

class Nucleic_codes
{
public:
    Nucleic_codes();
    QMap<char,char> getCharToByteMap();
    QMap<char, char> getByteToCharMap();
    char getCharFromByte(char byte);
    char getByteFromChar(char character);

    static const char dash = 1 << 4;
    static const char  A= 1 << 3;
    static const char  C= 1 << 2;
    static const char  G= 1 << 1;
    static const char T= 1 << 0;
    static const char all= dash | A | C | G | T;
private:
    QMap<char, char> charToByteMap;
    QMap<char, char> byteToCharMap;
    void fillCharToByteMap();
    void fillByteToCharMap();
};

#endif // NUCLEIC_CODES_H
