#ifndef PINDEL2VCF_H
#define PINDEL2VCF_H

#include <string>

void createComplement( const std::string& dna, std::string& complement );
char complementBase( char inputBase );

#endif

