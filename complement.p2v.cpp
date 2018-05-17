#include <sstream>
#include <vector>
#include <string> 
#include <iostream>
#include <fstream>
#include <map>
#include "GetPot"
#include "p2v.h"

using std::string;
using std::ios_base;
using std::ostream;
using std::cin;
using std::cout;
using std::cerr;
using std::ifstream;
using std::map;
using std::stringstream;
using std::vector;
using std::make_pair;
using std::endl;


typedef string Header;
typedef string Contig;

const int LINELENGTH = 60;

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


void print_revcomp(Header const& header, Contig const& seg, ostream& out = cout) {
    out << header << "\n";
    
    Contig comp(seg.rbegin(),seg.rend());
//  transform(comp.begin(),comp.end(), comp.begin(), complement);
		createComplement( seg, comp );
//	  out << comp << "\n";
    size_t i = 0;
    size_t stop = comp.length()/LINELENGTH + ((comp.length()%LINELENGTH)?1:0);
    
    while(i < stop)
       out << comp.substr(i++*LINELENGTH,LINELENGTH) << "\n";
}

int main (int argc, char **argv) {
  ios_base::sync_with_stdio(false);

  Contig sequence; 
  Header header;
	GetPot   cl(argc, argv);
	ifstream index_file, fasta_file;
	map< string, int > fa_map;
	string str_buf, line;
	int    int_buf;
	struct tokens {
		string ctg;
		int offset;
		int trash;
	};
	tokens buffy;
	//vector<string> tokens;

	const string fa  = cl("fasta", "default");
	const string ctg_id = cl( "ctg", "" );

	const string fai = fa + ".fai";

	index_file.open( fai, ifstream::in );

	if ( ! index_file.good() ) { 
		cerr << "failed to open: " << fai << "\n";
		exit(1);
	}
	else {
		while ( getline( index_file, line ) ) { 
			stringstream ss(line);
			ss >> buffy.ctg;
			ss >> buffy.trash;
			ss >> buffy.offset;
			fa_map.insert( make_pair( buffy.ctg, buffy.offset ) );

		}
	}

	map<string,int>::iterator it = fa_map.find( ctg_id );

	if ( it == fa_map.end() ) { 
		cerr << "not found: " << ctg_id << "\n";
		exit( 1 );
	}

	fasta_file.open( fa, ifstream::in );

	// if !c++11 or higher, need to fasta_file.clear() first in case eof was hit
	fasta_file.seekg( it->second );

  while (getline(fasta_file, line) && ! (line[0] == '>') ) {
		sequence += line;
  }
	if ( ! sequence.empty() ) { 
		header = ">" + ctg_id;
  	print_revcomp(header, sequence);
	}

  return 0;
}
