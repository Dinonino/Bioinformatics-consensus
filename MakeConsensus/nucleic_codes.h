#ifndef NUCLEIC_CODES_H
#define NUCLEIC_CODES_H
#include <map>

using namespace std;


class Nucleic_codes
{
public:
    Nucleic_codes();
    map<char,char> getCharToByteMap();
    map<char, char> getByteToCharMap();
    char getCharFromByte(char byte);
    char getByteFromChar(char character);

    static const char dash = 1 << 4;
    static const char A= 1 << 3;
    static const char C= 1 << 2;
    static const char G= 1 << 1;
    static const char T= 1 << 0;
private:
    map<char, char> charToByteMap;
    map<char, char> byteToCharMap;
    void fillCharToByteMap();
    void fillByteToCharMap();
};

#endif // NUCLEIC_CODES_H
