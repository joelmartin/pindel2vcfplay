#include <string> 
#include <iostream>
#include "p2v.h"

using std::string;
using std::ios_base;
using std::ostream;
using std::cin;
using std::cout;

typedef string Header;
typedef string Segment;

const int LINELENGTH = 60;

/* returns the complementary DNA-base of base 'inputbase' */
/*
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
*/

/* 'createComplement' creates the complement of a DNA string. */ // HEY!!! don't forget to complement CpG's properly!
/*
void createComplement( const string& dna, string& complement ) {
   int dnaLength = dna.length();
   complement = "";
   for (int position=dnaLength-1; position>=0; position-- ) {
      complement += complementBase( dna[ position ] );
   }
}
*/

void print_revcomp(Header const& header, Segment const& seg, ostream& out = cout) {
    out << header << "\n";
    
    Segment comp(seg.rbegin(),seg.rend());
    //transform(comp.begin(),comp.end(), comp.begin(), complement);
		createComplement( seg, comp );
//	  out << comp << "\n";
    size_t i = 0;
    size_t stop = comp.length()/LINELENGTH + ((comp.length()%LINELENGTH)?1:0);
    
    while(i < stop)
       out << comp.substr(i++*LINELENGTH,LINELENGTH) << "\n";
}


int main () {
  ios_base::sync_with_stdio(false);

  Segment line, segment; 
  Header header;

  while (getline(cin, line))
  {
      if (line[0] == '>')
      {
          if (! segment.empty())
          print_revcomp(header, segment);
          header = line;
          segment.clear();
      }
      else
          segment += line;
  }
  print_revcomp(header, segment);

  return 0;
}
