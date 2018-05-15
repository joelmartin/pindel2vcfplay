#include <string>
#include <iostream>
#include <fstream>
#include "p2v_chromosome.h"

using std::string;
using std::ifstream;
using std::cout;
using std::endl;
using std::getline;

bool g_normalBaseArray[256];

bool normalBase(char ch )
{
   return g_normalBaseArray[ ch ];
}

void makeStrangeBasesN(string& dna)
{
   int chromLength = dna.size();
   for (int position=0; position<chromLength; position++ ) {
      if (!normalBase(dna[ position ])) {
         dna[ position ] ='N';
      }
   }
}

void initBaseArray()
{
   for (int i=0; i<256; i++) {
      g_normalBaseArray[i] = false;
   }
   g_normalBaseArray['A'] = true;
   g_normalBaseArray['C'] = true;
   g_normalBaseArray['G'] = true;
   g_normalBaseArray['T'] = true;
   g_normalBaseArray['N'] = true;
}

/** 'convertToUppercase' returns the input string in full uppercase. */
string convertToUppercase( const string& inputString)
{
   string outputString = inputString;
   for (int i=0; i<outputString.length(); i++ ) {
      outputString[ i ] = toupper( outputString[ i ] );
   }
   return outputString;
}

void Chromosome::removeFromMemory() {
      cout << "Removing chromosome " << d_identifier << " from memory.\n";
      delete d_sequence;
      d_sequence=new string("");
}
/* 'Chromosome::readFromFile' reads in the reference chromosome sequence from the FASTA file. */
void Chromosome::readFromFile()
{
   ifstream referenceFile( d_fastaFilename.c_str() );
   referenceFile.clear(); // reset file for reading from the start again
   referenceFile.seekg(0);
   if (referenceFile.fail()) {
      cout << "Cannot open reference file. Exiting.\n";
      exit(EXIT_FAILURE);
   }
   string refLine, refName, currentLine;
   string tempChromosome = "";

   getline(referenceFile,refLine); // FASTA format always has a first line with the name of the reference in it
   // loop over each chromosome
   bool targetChromosomeRead = false;
   do {
      int counter=1;
      refName = "";
      do {
         refName += refLine[ counter++ ];
      } while ( counter<refLine.size() && (refLine[ counter ] != ' ') && (refLine[ counter ] != '\t') && (refLine[ counter ] != '\n') );

      if (refName == d_identifier ) {
         cout << "Reading chromosome " << refName << " into memory." << endl;
         targetChromosomeRead = true;
         tempChromosome+="N"; // for 1-shift
      }
      getline(referenceFile,currentLine);

      while (!referenceFile.eof() && currentLine[0]!='>') {
         if (refName == d_identifier ) {
            tempChromosome += convertToUppercase( currentLine );
         }
         getline(referenceFile,currentLine);
      }
      makeStrangeBasesN(tempChromosome);
      d_sequence = new string( tempChromosome );
      refLine = currentLine;
   } while (!referenceFile.eof() && !targetChromosomeRead);
}

const string* Chromosome::getChromPtr()
{
   if (d_sequence == NULL) { // sequence not read into memory
      readFromFile();
      //cout << "Have read " << d_sequence;
   }
   return d_sequence;
}


