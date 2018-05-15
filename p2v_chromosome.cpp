#include <string>
#include <iostream>
#include <fstream>
#include "p2v_chromosome.h"

using std::string;
using std::ifstream;
using std::cout;
using std::endl;

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
      char ch;
      referenceFile.get( ch );
      while (!referenceFile.eof() && ch != '>') {
         if (refName == d_identifier ) {
            char niceCh = toupper( ch );
            if (niceCh >='A' && niceCh<= 'Z' ) { // skip all spaces and tab stops etc.
               tempChromosome += niceCh;
            }
         }
         referenceFile.get( ch );
      }
      makeStrangeBasesN(tempChromosome);
      d_sequence = new string( tempChromosome );
      referenceFile.putback( ch );
      getline(referenceFile,refLine); // FASTA format always has a first line with the name of the reference in it
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


