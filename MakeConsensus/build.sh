mkdir bin

g++ -std=c++11 read.cpp unitig.cpp nucleic_codes.cpp  consensus.cpp realigner.cpp main.cpp -o bin/bio.out
