#include "nucleic_codes.h"


    const unsigned char Nucleic_codes::dash = 1 << 4;
    const unsigned char Nucleic_codes::A= 1 << 3;
    const unsigned char Nucleic_codes::C= 1 << 2;
    const unsigned char Nucleic_codes::G= 1 << 1;
    const unsigned char Nucleic_codes::T= 1 << 0;
    const unsigned char Nucleic_codes::all= dash | A | C | G | T;

Nucleic_codes::Nucleic_codes()
{
    fillCharToByteMap();
    fillByteToCharMap();
}

/**
 * @brief Nucleic_codes::fillCharToByteMap Metoda za popunjavanje mape koja sadrži podatke o sadržanim bajtovima na temelju danog znaka
 */
void Nucleic_codes::fillCharToByteMap(){
    charToByteMap['A']= A;
    charToByteMap['C']= C;
    charToByteMap['G']= G;
    charToByteMap['T']= T;
    charToByteMap['-']=dash;

    charToByteMap['R']= A | G;
    charToByteMap['Y']= C | T;
    charToByteMap['M']= C | A;
    charToByteMap['K']= T | G;
    charToByteMap['W']= T | A;
    charToByteMap['S']= C | G;
    charToByteMap['B']= C | T | G;
    charToByteMap['D']= A | T | G;
    charToByteMap['H']= A | T | C;
    charToByteMap['V']= A | C | G;
    charToByteMap['N']= A | C | G | T;

    charToByteMap['E']= A | dash;
    charToByteMap['F']= G | dash;
    charToByteMap['I']= C | dash;
    charToByteMap['J']= T | dash;
    charToByteMap['L']= A | G | dash;
    charToByteMap['O']= C | T | dash;
    charToByteMap['P']= C | A | dash;
    charToByteMap['U']= C | G | dash;
    charToByteMap['Z']= T | A | dash;
    charToByteMap['X']= T | G | dash;

    charToByteMap['!']= C | T | G | dash;
    charToByteMap['#']= A | T | G | dash;
    charToByteMap['$']= A | T | C | dash;
    charToByteMap['%']= A | C | G | dash;
    charToByteMap['&']= A | C | G | T | dash;

}

/**
 * @brief Nucleic_codes::fillByteToCharMap Metoda za popunjavanje mape koja sadrži podatke o sadržanim znaku na temelju danog bajta.
 */
void Nucleic_codes::fillByteToCharMap(){
    byteToCharMap[A]='A';
    byteToCharMap[C]='C';
    byteToCharMap[T]='T';
    byteToCharMap[G]='G';
    byteToCharMap[dash]='-';


    byteToCharMap[A | G]='R';
    byteToCharMap[C | T]='Y';
    byteToCharMap[C | A]='M';
    byteToCharMap[T | G]='K';
    byteToCharMap[T | A]='W';
    byteToCharMap[C | G]='S';

    byteToCharMap[C | T | G]='B';
    byteToCharMap[A | T | G]='D';
    byteToCharMap[A | T | C]='H';
    byteToCharMap[A | C | G]='V';
    byteToCharMap[A | C | G | T]='N';

    byteToCharMap[A | dash]='E';
    byteToCharMap[C | dash]='F';
    byteToCharMap[G | dash]='I';
    byteToCharMap[T | dash]='J';
    byteToCharMap[A | G | dash]='L';
    byteToCharMap[C | T | dash]='O';
    byteToCharMap[C | A | dash]='P';
    byteToCharMap[C | G | dash]='U';
    byteToCharMap[T | A | dash]='Z';
    byteToCharMap[T | G | dash]='X';

    byteToCharMap[C | T | G | dash]='!';
    byteToCharMap[A | T | G | dash]='#';
    byteToCharMap[A | T | C | dash]='$';
    byteToCharMap[A | C | G | dash]='%';
    byteToCharMap[A | C | G | T | dash]='&';

}

/**
 * @brief Nucleic_codes::getCharToByteMap Metoda za dohvaćanje mape koja sadrži podatke o sadržanim bajtovima na temelju danog znaka
 * @return Mapa koja sadrži podatke o sadržanim bajtovima na temelju danog znaka
 */
map<char,unsigned char> Nucleic_codes::getCharToByteMap(){
    return charToByteMap;
}

/**
 * @brief Nucleic_codes::getByteToCharMap Metoda za dohvaćanje mape koja sadrži podatke o sadržanim znaku na temelju danog bajta
 * @return Mape koja sadrži podatke o sadržanim znaku na temelju danog bajta
 */
map<unsigned char,char> Nucleic_codes::getByteToCharMap(){
    return byteToCharMap;
}

/**
 * @brief Nucleic_codes::getCharFromByte Metoda koja na temelju primljenog bajta vraća odgovarajući znak iz mape sa znakovima.
 * @param byte Bajt
 * @return Odgovarajući znak
 */
char Nucleic_codes::getCharFromByte(unsigned char byte){
    return byteToCharMap[byte];
}

/**
 * @brief Nucleic_codes::getByteFromChar Metoda koja na temelju primljenog znaka vraća njegovu odgovarajuću vrijednost u bajtovima.
 * @param character Znak
 * @return Odgovarajući bajt.
 */
unsigned char Nucleic_codes::getByteFromChar(char character){
    return charToByteMap[character];
}
