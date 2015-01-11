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

    static const char dash;
    static const char A;
    static const char C;
    static const char G;
    static const char T;
    static const char all;
private:
    map<char, char> charToByteMap;
    map<char, char> byteToCharMap;
    void fillCharToByteMap();
    void fillByteToCharMap();
};



#endif // NUCLEIC_CODES_H
