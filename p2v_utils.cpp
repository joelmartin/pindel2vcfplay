#include <string>
#include "p2v.h"

using std::string;

/* returns the complementary DNA-base of base 'inputbase' */
char complementBase( char inputBase ) {
	switch (inputBase) {
	case 'A':
		return 'T';
	case 'C':
		return 'G';
	case 'G':
		return 'C';
	case 'T':
		return 'A';
	default:
		return 'N';
	}
}

/* 'createComplement' creates the complement of a DNA string. */ // HEY!!! don't forget to complement CpG's properly!
void createComplement( const string& dna, string& complement ) {
	int dnaLength = dna.length();
	complement = "";
	for (int position=dnaLength-1; position>=0; position-- ) {
		complement += complementBase( dna[ position ] );
	}
}
