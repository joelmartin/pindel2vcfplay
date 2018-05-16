/* vcfcreator.cpp

 Transforms indel and SV-files created by Pindel into VCF4.0 format according to 1000-genome SV specifications.
 This version sorts the records before outputting; may need to make version that uses external memory later!

 Created by Eric-Wubbo Lameijer, section of Molecular Epidemiology, Leiden University Medical Center, March 3rd, 2011.
 e.m.w.lameijer@gmail.com
 +31(0)71-5 125 831

 Version 0.6.3 [February 19th, 2014] Clearer text on usage of -P option
 Version 0.6.2 [December 12th, 2014] Now robust against fasta files that have non-standard line lengths (C++'s getline does not work well on lines of over a million characters)
 Version 0.6.1 [December 12th, 2014] Now has special code to recognize lines that contain SV-data, instead of relying on indirect establishment of their identity from context
		(for example: seeing that the previous line starts with '#', which fails for a grepped Pindel output file.)
 Version 0.6.0 [April 21nd, 2014] automatic warning when finding oversized SVs - and set -co by default as a safeguard. Thanks to user Hong Ching Lee for reporting this problem.
 Version 0.5.9 -ho -he for cutoffs
 Version 0.5.8 [April 28th, 2013] Added "-co" option for compact output; -co followed by an integer indicates
 the longest variation displayed with full base sequence. If, for example, the user gives as command line parameter
 "-co 10", then if either the reference or the alternative allele exceeds 10 bases, the output will be transformed
 into the format chrom pos firstrefbase <SVType>.
 Version 0.5.7 [February 1st, 2013] Now updated my contact information.
 Version 0.5.6 [February 1st, 2013]  Quickfix since user detected empty labels. Note: since no input was delivered that enabled me to reproduce the bug, I could not find the true cause. So problems may still crop up in the future.
 Version 0.5.5 [December 10th, 2012] Modified the code so that LI / -P now gives correct GT:AD instead of GT:RD:AD output. Also fixed small bug in creating ALT of NT-inversions
 Version 0.5.4 [December 6th, 2012] Error found by David Hannah on using the -P option debugged
 Version 0.5.3 [November 8th, 2012] Now should genotype newer pindel output as 0/0, 1/1 or 0/1 based on apparent balance between alleles
 Version 0.5.2 [November 8th, 2012] Bugs removed from 0.5.1 and added GT:RD:AD output to show reference coverage
 Version 0.5.1 [October 31st, 2012] Now also compatible with Pindel's new output format
 Version 0.5.0 [July 9th, 2012] removed small error that caused the -c option to process all other chromosomes as well...
 Version 0.4.9 [July 2nd, 2012] added -P option to allow fusing all output files of one pindel run
 Version 0.4.8 [July 2nd, 2012] debugged .:10 genotype (due to 'dual-encoding' of genotype, while Pindel not having actual chromosome data
 Version 0.4.7 [June 21th, 2012] LI now also have proper labels <as per updated Pindel>
 Version 0.4.6 [June 20th, 2012] END-position now proper according to http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41, access date June 20, 2012
 Version 0.4.5 [June 20th, 2012] displays proper warnings when using the -G option
 Version 0.4.4 [June 19th, 2012] now also prints inversions correctly in GATK-format (when using the -G option)
 Version 0.4.3 [June 19th, 2012] adding equilength replacement calls to -G output (2)
 Version 0.4.2 [June 19th, 2012] adding equilength replacement calls to -G output (1)
 Version 0.4.1 [June 18th, 2012] Now automatically sets the end of long insertions to occur after the beginning...
 Version 0.4.0 [June 18th, 2012] Added -G option to make output GATK-compatible
 Version 0.3.9 [June 8th, 2012] Added maximum total coverage to protect from duplicated regions
 Veriosn 0.3.8 [June 8th, 2012] Long insertions support debugged.
 Version 0.3.7 [June 8th, 2012] Next stage in debugging repeats/postindel options
 Version 0.3.6 [June 8th, 2012] Debugged int repeats/postindel options[1]
 Version 0.3.5 [June 1st, 2012] Refined microhomology/microsattelite sequencing to distinguish internal repeats and 'postindel' repeats
 Version 0.3.4 [May 21st, 2012] Alternative homology calling; now checks whether inserted/deleted sequence is part of a longer repetitive sequence
 Version 0.3.3 [May 2nd, 2012] Debugged -sb/-ss option: now works if 1 sample is selected
 Version 0.3.2 [May 1st, 2012] Adds -sb and -ss options to allow users to only count samples with sufficient individual support.
 Version 0.3.1 [March 27th, 2012] -c option now works!
 Version 0.3.0 [March 27th, 2012] now also works correctly for small inversions (as NT is then only "TG", gave wrong code, ACA ATGTGTG instead of ACA ATG
 Version 0.2.9 [March 27th, 2012] Adds an option for counting samples only if they are balanced with at least N reads on each side
 Version 0.2.8 [March 2nd, 2012] Now does not skip chromosome when compiling long insertions
 Version 0.2.7 [March 2nd, 2012] Skip chromosomes in which Pindel has not called any SVs
 Version 0.2.6 [March 2nd, 2012] Fixed problem in which chromosome was not properly freed from memory
 Version 0.2.5 [March 2nd, 2012] Fixed segmentation fault error when there are no reads in a chromosome
 Version 0.2.4 [February 6th, 2012] Can save memory at the cost of time by implementing a window-size option.
 Version 0.2.3 [February 6th, 2012] More memory-efficient by not explicitly storing REF and ALT strings, but by using recalculation instead of storage
 Version 0.2.1 [January 3rd, 2012] Now also correctly reads in fasta files which do not only contain uppercase characters.
 Version 0.2.0 [November 10th, 2011. The support of an indel, previously called "DP" is now more appropriately called "AD".
 Also, genotypes are somewhat more correctly indicated with . and 1/.
 Also, replaces -1 by . to indicate unknown number of fields in the declarations in the header
 Version 0.1.9 [August 19th, 2011. To save memory, now reads in chromosomes only when needed, so doesn't put the entire genome in memory at once.
 */

/* CONTENTS

 PREFACE: include files and global constants
 CHAPTER 1. General utilities for DNA-string manipulation (reverse-complementing a string)
 CHAPTER 2. Defining the parameters and the 'Parameter class' to handle them.
 */


/*** PREFACE: include files and global constants ->***/

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <vector>

#include <ctype.h> // tolower
#include <math.h> //sqrt
#include <stdlib.h> // for atoi
#include <time.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>

#include "p2v.h"
#include "p2v_input_reader.h"
#include "p2v_parameters.h"
#include "p2v_chromosome.h"
#include "p2v_genome.h"
#include "p2v_svdata.h"

const int FIRST_SAMPLE_INDEX = 32; // index of first sample name

//jam.note to remove using namespace std with exlicit and limited usings of std
using namespace std;
using std::getline;
using std::cout;

string g_versionString = "0.6.3";
string g_programName = "pindel2vcf";


int  g_sizeToWarnFor = 1000000;
bool pindel024uOrLater = false;

int getsizeToWarnFor() {
	return( g_sizeToWarnFor);
}
bool getpindel024uOrLater () {
	return( pindel024uOrLater );
}

//bool pindel024uOrLater = false;

/* the global parameter g_par stores the values of all parameters set by the user; the Parameter class, in contrast, is for user-friendly IO. */
ParameterSettings g_par;

bool getGATKCompatible( ) {
	return g_par.gatkCompatible;
}

int getCompactOutput () {
	return g_par.compactOutput;
}

bool getSomatic () {
	return g_par.somatic;
}

int getMinCoverage () {
	return g_par.MinCoverage;
}

double getHetCutoff () {
	return g_par.HetCutoff;
}

double getHomCutoff () {
	return g_par.HomCutoff;
}

/*** END OF PREFACE: include files and global constants ***/

vector<Parameter*> parameters;

/* 'createHeader' writes a VCF4.0 (1k genome SV) compatible header to output VCF-file "outFile". */
void createHeader(ofstream &outFile, const string& sourceProgram, const string& reference, const set<string>& samples) {
	// file format
	outFile << "##fileformat=VCFv4.0\n";

	// date of file creation
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	outFile << "##fileDate=" << g_par.referenceDate << endl;

	// source
	outFile << "##source=" << sourceProgram << endl;

	// reference
	outFile << "##reference=" << reference << endl;

	// info fields (selected from http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/VCF%20%28Variant%20Call%20Format%29%20version%204.0/encoding-structural-variants)
	outFile << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" << endl;
	outFile << "##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">" << endl;
	outFile << "##INFO=<ID=PF,Number=1,Type=Integer,Description=\"The number of samples carry the variant\">" << endl;
	outFile << "##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">" << endl;
	outFile << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << endl;
	outFile << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << endl;
	outFile << "##INFO=<ID=NTLEN,Number=.,Type=Integer,Description=\"Number of bases inserted in place of deleted code\">" << endl;
	outFile << "##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">" << endl;
	//outFile << "##ALT=<ID=DEL,Description=\"Deletion\">" << endl; /*EWL040311: probably not needed, as our calls are precise so we should rather give the exact replacing sequence instead of a label. */
	//outFile << "##ALT=<ID=DUP,Description=\"Duplication\">" << endl;
	//outFile << "##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">" << endl;
	//outFile << "##ALT=<ID=INS,Description=\"Insertion of novel sequence\">" << endl;
	//outFile << "##ALT=<ID=INV,Description=\"Inversion\">" << endl;
	//outFile << "##ALT=<ID=CNV,Description=\"Copy number variable region\">" << endl;
	outFile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
	if (pindel024uOrLater) {
		outFile << "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Reference depth, how many reads support the reference\">" << endl;
	}
	outFile << "##FORMAT=<ID=AD,Number=2,Type=Integer,Description=\"Allele depth, how many reads support this allele\">" << endl;

	// headers of columns
	outFile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
	if (samples.size()>0) {
		outFile << "\tFORMAT";
		for (set<string>::iterator counter=samples.begin(); counter!=samples.end(); counter++ ) {
			outFile << "\t" << *counter;
		}
	}
	outFile << "\n";
}

/* 'Pair' contains a pair of strings, outputted as "str1:str2". e*/
class Pair {
	friend ostream& operator<<(ostream& os, const Pair& pair);

public:
	Pair(const string& first, const string& second) {
		d_first=first;
		d_second=second;
	}

private:
	string d_first, d_second;

};

ostream& operator<<(ostream& os, const Pair& pair) {
	os << pair.d_first << "," << pair.d_second;
    return os; //hack for xcode compiler.jam
}

// svdata was starting here

/* 'fetchElement' returns the "index"th element of "infile"; so the first element if index is 1, the second if index is 2, etc. */
string fetchElement( istream& instream, const int index ) {
	string element = "";
	for ( int currentCounter=0; currentCounter<index; currentCounter++ ) {
		instream >> element;
	}
	return element;
}

/* 'countElements' counts the number of elements in stream "instream" */
int countElements( istream& instream ) {
	string element = "";
	int counter = 0;
	while (! instream.fail() ) {
		instream >> element;
		counter++;
	}
	return counter;
}

void showSet( set<string> aSet ) {

	set<string>::iterator index;
	int counter=1;
	for (index=aSet.begin(); index!=aSet.end(); index++ ) {
		cout << counter++ << ". " << *index << endl;
	}
}

/** 'isPindelSVIdentifier' returns whether the string passed to it represents a valid SV type (insertion, deletion, etc.)
	 NOTE: this is the SV identifier as used in pindel output files, pindel2vcf uses other SV types internally,
	 like 'INS' or 'RPL' */
bool isPindelSVIdentifier( string identifier ) {
	if (identifier == "D" || identifier == "I" || identifier == "LI" || identifier == "TD" || identifier == "INV" ) {
		return true;
	}
	else {
		return false;
	}
}

/** 'isSVSummarizingLine' checks whether 'line' is a SV summarizing line and therefore one that pindel2vcf needs to process
	 (other lines can be skipped) */
bool isSVSummarizingLine( string line ) {
	// first check if the line is a true summarizing name by counting if it has sufficient elements
	// (this is also useful to prevent reading of the SV type at the second position of the line to
	// go wrong if there is no second element)
	stringstream streamForCounting;
	streamForCounting << line;
	int elementsInLine = countElements( streamForCounting );
	if (elementsInLine < 2 ) {
		return false;
	}

	// there are at least two elements. Check if the second one indicates a SV-type
	stringstream lineStream;
	lineStream << line;
	string item = fetchElement( lineStream, 2 );
	return isPindelSVIdentifier( item );
}

void getSampleNamesAndChromosomeNames(InputReader& pindelInput, set<string>& sampleNames, set<string>&chromosomeNames) {
	//cout << "DEBUG:start GetSampleNamesAndChromosomeNames\n";
	string line;
	unsigned int counter=0;

	// loop over the input file, registering each 'interesting' pindel output line (pindel output has the structure (ABCDn)m
	// A is a ##### separator line, B the interesting line (with sample names and coverage), C the line with the reference DNA
	// sequence (starting with ACTG) and the D-lines start with spaces and give the individual reads that support an event
	// HOWEVER, it is also possible that an user greps only the 'useful' output lines.
	while (!pindelInput.eof()) {

		// find the next 'B-line' (the line containing the sample names, coverage etc.)
		do {
			line = pindelInput.getLine();
		}
		while (!pindelInput.eof() && !isSVSummarizingLine( line ));

		if (pindelInput.eof()) {
			//cout << "DEBUG:end GetSampleNamesAndChromosomeNames\n";
			return;
		}

		// 'line' should contain a B-line now. Analyze it.
		counter++;
		stringstream lineStream;
		lineStream << line;
		stringstream streamForCounting;
		streamForCounting << line;
		int elementsInLine = countElements( streamForCounting );
		string svType = fetchElement( lineStream, 2 );


		// 'LI' types don't have the same format as normal Pindel lines.
		if ( svType.compare("LI")==0 ) {
			string chromosomeName = fetchElement( lineStream, 2 );
			chromosomeNames.insert( chromosomeName );
			string firstSampleName = fetchElement( lineStream, 7);
			sampleNames.insert( firstSampleName );
			string newSampleName = fetchElement( lineStream, 5 );
			while (!lineStream.fail()) {
				sampleNames.insert( newSampleName );
				newSampleName = fetchElement( lineStream, 5 );
			}
			continue;
		}
		string chromosomeName = fetchElement( lineStream, 6 );
		chromosomeNames.insert( chromosomeName );
		//		cout << "Studying chromosome " << chromosome << endl;

		// 8 = 2+6, so corrects for previous reads
		int numberOfSamples = atoi( fetchElement( lineStream, FIRST_SAMPLE_INDEX - 12 ).c_str() );
		string firstSampleName = fetchElement( lineStream, 4 );
		if ( firstSampleName != "" ) {
			sampleNames.insert( firstSampleName );
		}

		//cout << "ElInLine: " << elementsInLine << ", FSINDEX: " << FIRST_SAMPLE_INDEX << ", NoS=" << numberOfSamples << endl;
		if (elementsInLine> FIRST_SAMPLE_INDEX + 5* numberOfSamples ) {
			pindel024uOrLater = true;
		}
		/*      else { // *** This code seems to give trouble with some pindel output; pindel output format not consistent?
		 pindel024uOrLater = false;
		 }*/
		int numberOfElementsPerSample = ( pindel024uOrLater ? 7 : 5 );
		string newSampleName = fetchElement( lineStream, numberOfElementsPerSample );
		while (!lineStream.fail()) {
			if ( newSampleName != "" ) {
				sampleNames.insert( newSampleName );
			}
			newSampleName = fetchElement( lineStream, numberOfElementsPerSample );
		}
	}

	//cout << "DEBUG:end GetSampleNamesAndChromosomeNames\n";
}

template<class T> void showVector( vector<T> vect ) {
	cout << "Printing vector: " << endl;
	for (int index=0; index<vect.size(); index++ ) {
		cout << index << " " << vect[ index ] << endl;
	}
	cout << "End of printing.\n";
}

/* 'convertIndelToSVdata' converts insertions and deletions to nicely formatted SV-data. */
void convertIndelToSVdata( InputReader& pindelInput, map< string, int>& sampleMap, Genome& genome, SVData& svd, const string& targetChromosomeID) {
	string line;
	svd.setGenome( genome );
	do {
		line = pindelInput.getLine();
	}
	while (!pindelInput.eof() && !isSVSummarizingLine( line ));

	if (pindelInput.eof()) {
		return;
	}

	stringstream lineStream;
	lineStream << line;
	string svType = fetchElement( lineStream, 2 ); // to 2
	if ( svType.compare("LI") == 0 ) {
		svd.setSVtype("INS");
		svd.setSVlen( 0 );
		string chromosomeID = fetchElement( lineStream, 2);
		const string* reference = genome.getChromosome( chromosomeID );
		if ( reference== NULL ) {
			cout << "Error! Reference chromosome \"" << chromosomeID << "\" not found!" << endl;
			exit(EXIT_FAILURE);
		}
		svd.setChromosome( chromosomeID );
		if (chromosomeID!=targetChromosomeID) {
			return;
		}
		int beforeStartPos = atoi( fetchElement( lineStream, 1 ).c_str() );
		svd.setPosition( beforeStartPos );
		int totalPlusSupport = atoi( fetchElement( lineStream, 2 ).c_str());
		int rightmostEndPos = atoi (fetchElement( lineStream, 1 ).c_str()); // now at position 14
		//cout << "plusSupport, righmostEndPos is " << plusSupport << ", " << rightmostEndPos << endl;
		svd.setEnd( rightmostEndPos );
		svd.setBPrange( beforeStartPos, rightmostEndPos );
		int totalMinSupport = atoi( fetchElement( lineStream, 2 ).c_str());
		// if the file has been created by a recent version of pindel, read in the extra elements
		string sampleName = fetchElement( lineStream, 1);
		int refSupportAtStartOfEvent = 0;
		int refSupportAtEndOfEvent = 0;

		/*if ( pindel024uOrLater ) {
		   refSupportAtStartOfEvent = atoi( fetchElement( lineStream, 1 ).c_str());
		   refSupportAtEndOfEvent = atoi( fetchElement( lineStream, 1 ).c_str());
		   }*/
		int totalRefSupport = max(refSupportAtStartOfEvent, refSupportAtEndOfEvent);
		int numberItemsUntilNextSupport = ( pindel024uOrLater ? 2 : 2 );
		int samplePlusSupport = atoi( fetchElement( lineStream, numberItemsUntilNextSupport ).c_str());
		int sampleMinSupport = atoi( fetchElement( lineStream, numberItemsUntilNextSupport ).c_str()); // now at position 35, total +supports sample 1
		//int count=0;
		while (!lineStream.fail()) {
			if (sampleMap.find( sampleName )==sampleMap.end() ) {
				cout << "Error: could not find sample " << sampleName << endl;
			}
			else {
				int sampleID = sampleMap[ sampleName ];
				svd.addGenotype( sampleID, samplePlusSupport , sampleMinSupport, totalRefSupport );
			}
			/*if ( pindel024uOrLater ) {
			 refSupportAtStartOfEvent = atoi( fetchElement( lineStream, 1 ).c_str() );
			 refSupportAtEndOfEvent = atoi( fetchElement( lineStream, 1 ).c_str());
			 }*/
			int totalRefSupport = max(refSupportAtStartOfEvent, refSupportAtEndOfEvent);
			sampleName = fetchElement( lineStream, 1); // for unique support, 2->1,
			samplePlusSupport = atoi( fetchElement( lineStream, numberItemsUntilNextSupport ).c_str()); // for unique support, 2->1
			sampleMinSupport = atoi( fetchElement( lineStream, numberItemsUntilNextSupport ).c_str()); // now at position 33, total +supports sample 1
		}
		return;
	}
	svd.setSVlen( fetchElement( lineStream, 1 ) ); // to 3

	// get number(s) of NT bases added (two numbers for inversions!)
	string numNTaddedStr = fetchElement( lineStream, 2 ); // to 5

	int numNTadded = atoi( numNTaddedStr.c_str() ); // should get first number
	bool simpleInversion = false;
	int numNTinvAdded=-1;
	//cout << "Printing " << numNTaddedStr << endl;
	if ( svType.compare("INV") == 0 ) { // two numbers separated by : instead of one
		//cout << "Found ':' at position " << numNTaddedStr.find(":") << endl;
		if (numNTaddedStr.find(":")==string::npos) {
			//cout << "Found simple inversion!\n";
			simpleInversion = true;
		}
		else {
			int separatorPos = numNTaddedStr.find(":");
			string secondNumber = numNTaddedStr.substr(separatorPos+1);
			numNTinvAdded = atoi( secondNumber.c_str() );
		}
	}

	string ntAdded = fetchElement( lineStream, 1 ); // to 6
	string ntInvAdded = "";

	// basically, there are two type of inversions:
	//	a) The 'alternatively called small deletions' INV 2 NT 2 "TG"
	// b) the regular inversions INV98 NT 0:60 "":"GCT"
	if ( svType.compare("INV") == 0 ) {
		if (ntAdded.find(":")==string::npos ) {
			simpleInversion = true;
		}
		else {
			int separatorPos = ntAdded.find(":");
			ntInvAdded = ntAdded.substr( separatorPos+2, numNTinvAdded ); // erases ""
			svd.setSecondNT( ntInvAdded );
			ntAdded = ntAdded.substr(0,separatorPos);
		}
	}
	ntAdded.erase(0,1); // erases opening "
	ntAdded.erase(numNTadded); // erases closing "
	if (!simpleInversion) {
		svd.setNT( ntAdded );
	}

	string chromosomeID = fetchElement( lineStream, 2); // now at position 8
	if (chromosomeID!=targetChromosomeID) {
		return;
	}
	const string* reference = genome.getChromosome( chromosomeID );
	//cout << "reference is " << *reference << endl;
	if ( reference== NULL ) {
		cout << "Error! Reference chromosome \"" << chromosomeID << "\" not found!" << endl;
		exit(EXIT_FAILURE);
	}
	svd.setChromosome( chromosomeID );
	int beforeStartPos = atoi( fetchElement( lineStream, 2 ).c_str() ); // pos 10
	svd.setPosition( beforeStartPos ); // now at position 10
	int leftmostEndPos = atoi( fetchElement( lineStream, 1 ).c_str()); // now at position 11
	int leftmostStartPos = atoi (fetchElement( lineStream, 2 ).c_str());  // at position 13
	int rightmostEndPos = atoi (fetchElement( lineStream, 1 ).c_str()); // now at position 14
	svd.setBPrange( leftmostStartPos, rightmostEndPos );
	svd.setEnd( leftmostEndPos );
	svd.setHomlen( rightmostEndPos - leftmostEndPos );
	string homSeq="";
	for (int position=leftmostEndPos; position<rightmostEndPos; position++ ) {
		homSeq += (*reference)[ position ];
	}
	svd.setHomseq( homSeq );
	if ( svType.compare("D")==0 ) {
		if (numNTadded==0 ) {
			svd.setSVtype( "DEL" );
			svd.setReplace( 0 );
		}
		else {   // some NT-bases added
			svd.setSVtype( "RPL" );
			svd.setReplace( numNTadded );
		}
	}
	else if ( svType.compare("I")==0 ) {
		svd.setSVtype("INS");
		svd.setReplace( 0 );
	}
	else if ( svType.compare("TD")==0 ) {
		svd.setSVtype("DUP:TANDEM");
		svd.setReplace( numNTadded );
	}
	else if ( svType.compare("INV") == 0 ) {
		svd.setSVtype("INV");
		if (simpleInversion) {
			svd.setReplace( 0, 0 );
		}
		else {
			svd.setReplace( numNTadded, numNTinvAdded );
		}
	}
	string sampleName = fetchElement( lineStream, 18);
	int refSupportAtStartOfEvent = 0;
	int refSupportAtEndOfEvent = 0;

	if ( pindel024uOrLater ) {
		refSupportAtStartOfEvent = atoi( fetchElement( lineStream, 1 ).c_str() );
		refSupportAtEndOfEvent = atoi( fetchElement( lineStream, 1 ).c_str() );
		//std::cout << "outside refSupportAtStartOfEvent: " <<  sampleName << "\t" << refSupportAtStartOfEvent << "\t" << "refSupportAtEndOfEvent: " << refSupportAtEndOfEvent << std::endl;
	}
	int totalRefSupport = max(refSupportAtStartOfEvent, refSupportAtEndOfEvent);
	//std::cout << "outside totalRefSupport: " << totalRefSupport << std::endl;
	int numberOfItemsUntilNextSupport = ( pindel024uOrLater ? 2 : 2 );
	int plusSupport = atoi( fetchElement( lineStream, numberOfItemsUntilNextSupport - 1 ).c_str()); // now at position 33, total +supports sample 1; for unique support 1->2
	int minSupport = atoi( fetchElement( lineStream, numberOfItemsUntilNextSupport ).c_str()); // now at position 35, total +supports sample 1
	int count=0;
	while (!lineStream.fail()) {
		if (sampleMap.find( sampleName )==sampleMap.end() ) {
			cout << "Error: could not find sample " << sampleName << endl;
		}
		else {
			int sampleID = sampleMap[ sampleName ];
			//std::cout << "Adding " << sampleID << "\t" << plusSupport << "\t" << minSupport << "\t" << totalRefSupport << std::endl;
			svd.addGenotype( sampleID, plusSupport , minSupport, totalRefSupport );
		}
		sampleName = fetchElement( lineStream, 2); // for unique support, 2->1,
		if ( pindel024uOrLater ) {
			refSupportAtStartOfEvent = atoi( fetchElement( lineStream, 1 ).c_str() );
			refSupportAtEndOfEvent = atoi( fetchElement( lineStream, 1 ).c_str() );
			//std::cout << "inside refSupportAtStartOfEvent: " << sampleName << "\t" << refSupportAtStartOfEvent << "\t" << "refSupportAtEndOfEvent: " << refSupportAtEndOfEvent << std::endl;
		}
		totalRefSupport = max(refSupportAtStartOfEvent, refSupportAtEndOfEvent);
		//std::cout << "insert totalRefSupport: " << totalRefSupport << std::endl;
		plusSupport = atoi( fetchElement( lineStream, numberOfItemsUntilNextSupport - 1 ).c_str()); // for unique support, 2->1
		minSupport = atoi( fetchElement( lineStream, numberOfItemsUntilNextSupport ).c_str()); // now at position 33, total +supports sample 1
	}
}

/* 'readReference' reads in the reference. */
void readReference( const string& referenceName, Genome& genome ) {
	ifstream referenceFile( referenceName.c_str() );
	if (referenceFile.fail()) {
		cout << "Cannot open reference file. Exiting.\n";
		exit(EXIT_FAILURE);
	}
	//reference="N"; // trick to make the reference automatically 1-positioned
	string refLine, refName, currentLine;

	getline(referenceFile,refLine); // FASTQ format always has a first line with the name of the reference in it
	// loop over each chromosome
	do {

		int counter=1;
		refName = "";
		do {
			refName += refLine[ counter++ ];
		}
		while ( counter<refLine.size() && (refLine[ counter ] != ' ') && (refLine[ counter ] != '\t') && (refLine[ counter ] != '\n') && (refLine[ counter ] != '\r'));
		cout << "Scanning chromosome: " << refName << endl;
		Chromosome newChrom( refName, referenceName );
		getline(referenceFile,currentLine);
		while (!referenceFile.eof() && currentLine[0]!='>') {
			getline(referenceFile,currentLine);
		}
		genome.addChromosome( newChrom );
		refLine = currentLine;
	}
	while (!referenceFile.eof());
	//cout << "DEBUG:Exiting reference scanning.\n";
}


/* 'createParameters' creates the default parameters that the VCF converter uses. */
void createParameters() {
	parameters.push_back(
	  new StringParameter( &g_par.reference, "-r", "--reference", "The name of the file containing the reference genome", true, "" ) );
	parameters.push_back(
	  new StringParameter( &g_par.referenceName, "-R", "--reference_name", "The name and version of the reference genome", true, "" ) );
	parameters.push_back(
	  new StringParameter( &g_par.referenceDate, "-d", "--reference_date", "The date of the version of the reference genome used", true, "" ) );
	parameters.push_back(
	  new StringParameter( &g_par.pindelfile, "-p", "--pindel_output", "The name of the pindel output file containing the SVs", false, "" ) );
	parameters.push_back(
	  new StringParameter( &g_par.pindelroot, "-P", "--pindel_output_root", "The root-name of the pindel output file; this will result in\n"
	                       "one big output file containing deletions, short and long insertions, tandem duplications and inversions.\n"
	                       "For example, if the pindel output files are called sample1_D, sample1_SI, sample1_TD etc. then -P sample1 would combine the\n"
	                       "information in all those sample files into one big vcf file.", false, "" ) );
	parameters.push_back(
	  new StringParameter( &g_par.vcffile, "-v", "--vcf", "The name of the output vcf-file (default: name of pindel output file +\".vcf\"", false, "" ) );
	parameters.push_back(
	  new StringParameter( &g_par.chromosome, "-c", "--chromosome", "The name of the chromosome (default: SVs on all chromosomes are processed)", false, "" ) );
	parameters.push_back(
	  new IntParameter( &g_par.windowSize, "-w", "--window_size", "Memory saving option: the size of the genomic region in a chromosome of which structural variants are calculated separately, in millions of bases (default 300, for memory saving 100 or 50 recommended)", false, 300 ) );
	parameters.push_back(
	  new IntParameter( &g_par.MinCoverage, "-mc", "--min_coverage", "The minimum number of reads to provide a genotype (default 10)", false, 10 ) );
	parameters.push_back(
	  new FloatParameter( &g_par.HetCutoff, "-he", "--het_cutoff", "The propertion of reads to call het (default 0.2)", false, 0.2 ) );
	parameters.push_back(
	  new FloatParameter( &g_par.HomCutoff, "-ho", "--hom_cutoff", "The propertion of reads to call het (default 0.8)", false, 0.8 ) );

	parameters.push_back(
	  new IntParameter( &g_par.minsize, "-is", "--min_size", "The minimum size of events to be reported (default 1)", false, 1 ) );
	parameters.push_back(
	  new IntParameter( &g_par.maxsize, "-as", "--max_size", "The maximum size of events to be reported (default infinite)", false, -1 ) );
	parameters.push_back(
	  new BoolParameter( &g_par.bothstrands, "-b", "--both_strands_supported", "Only report events that are detected on both strands (default false)", false, false ) );
	parameters.push_back(
	  new IntParameter( &g_par.minsuppSamples, "-m", "--min_supporting_samples", "The minimum number of samples an event needs to occur in with sufficient support to be reported (default 0)", false, 1 ) );
	parameters.push_back(
	  new IntParameter( &g_par.minsuppReads, "-e", "--min_supporting_reads", "The minimum number of supporting reads required for an event to be reported (default 1)", false, 1 ) );
	parameters.push_back(
	  new IntParameter( &g_par.maxSuppReads, "-f", "--max_supporting_reads", "The maximum number of supporting reads allowed for an event to be reported, allows protection against miscalls in due to segmental duplications or poorly mapped regions (default infinite)", false, -1 ) );
	parameters.push_back(
	  new IntParameter( &g_par.regionStart, "-sr", "--region_start", "The start of the region of which events are to be reported (default 0)", false, 0 ) );
	parameters.push_back(
	  new IntParameter( &g_par.regionEnd, "-er", "--region_end", "The end of the region of which events are to be reported (default infinite)", false, -1 ) );
	parameters.push_back(
	  new IntParameter( &g_par.maxInterRepeatNo, "-ir", "--max_internal_repeats", "Filters out all indels where the inserted/deleted sequence is a homopolymer/microsatellite of more than X repetitions (default infinite). For example: T->TCACACA has CACACA as insertion, which is a microsattelite of 3 repeats; this would be filtered out by setting -ir to 2", false, -1 ) );
	parameters.push_back(
	  new IntParameter( &g_par.compactOutput, "-co", "--compact_output_limit", "Puts all structural variations of which either the ref allele or the alt allele exceeds the specified size (say 10 in '-co 10') in the format 'chrom pos first_base <SVType>'", false, g_sizeToWarnFor ) );
	parameters.push_back(
	  new IntParameter( &g_par.maxInterRepeatLength, "-il", "--max_internal_repeatlength", "Filters out all indels where the inserted/deleted sequence is a homopolymers/microsatellite with an unit size of more than Y, combine with the option -ir. Default value of -il is infinite. For example: T->TCAGCAG has CAGCAG as insertion, which has the fundamental repetitive unit CAG of length 3. This would be filtered out if -il has been set to 3 or above, but would be deemed 'sufficiently unrepetitive' if -il is 2", false, -1 ) );
	parameters.push_back(
	  new IntParameter( &g_par.maxPostRepeatNo, "-pr", "--max_postindel_repeats", "Filters out all indels where the inserted/deleted sequence is followed by a repetition (of over X times) of the fundamental repeat unit of the inserted/deleted sequence. For example, T->TCACA would usually be a normal insertion, which is not filtered out, but if the real sequence change is TCACACA->TCACACACACA, it will be filtered out by -pr of 1 or above, as the fundamental repeat unit of the inserted sequence (CA) is repeated more than one time in the postindel sequence [indel sequence CACA, postindel sequence CACACA]. Note: when CAC is inserted next to ACACAC, the repeat sequence is recognized as CA, even though the 'postrepeat' sequence is ACACAC", false, -1 ) );
	parameters.push_back(
	  new IntParameter( &g_par.maxPostRepeatLength, "-pl", "--max_postindel_repeatlength", "Filters out all indels where the inserted/deleted sequence is followed by a repetition of  the fundamental repeat unit of the inserted/deleted sequence; the maximum size of that 'fundamental unit' given by the value of -pl (default infinite) For example: TCAG->TCAGCAG has insertion CAG and post-insertion sequence CAG. This insertion would be filtered out if -pl has been set to 3 or above, but would be deemed 'sufficiently unrepetitive' if -pl is 2", false, -1 ) );
	parameters.push_back(
	  new BoolParameter( &g_par.onlyBalancedSamples, "-sb", "--only_balanced_samples", "Only count a sample as supporting an event if it is supported by reads on both strands, minimum reads per strand given by the -ss parameter. (default false)", false, 0 ) );
//   parameters.push_back(
	//    new BoolParameter( &g_par.somatic, "-so", "--somatic_p", "compute somatic p value when two samples are present, assume the order is normal and tumor. (default false)", false, 0 ) );

	parameters.push_back(
	  new IntParameter( &g_par.minimumStrandSupport, "-ss", "--minimum_strand_support", "Only count a sample as supporting an event if at least one of its strands is supported by X reads (default 1)", false, 1 ) );
	parameters.push_back(
	  new BoolParameter( &g_par.gatkCompatible, "-G", "--gatk_compatible", "calls genotypes which could either be homozygous or heterozygous not as ./1 but as 0/1, to ensure compatibility with GATK", false, false ) );
	parameters.push_back(
	  new BoolParameter( &g_par.showHelp, "-h", "--help", "Print the help of this converter", false, false ) );
}

/* 'findParameter' returns the index of the parameter with name 'name'; -1 if not found.*/
int findParameter(string name) {
	for (int parameterCounter=0; parameterCounter<parameters.size(); parameterCounter++ ) {
		if (parameters[ parameterCounter ]->hasName( name ) ) {
			return parameterCounter;
		}
	}
	return -1;
}

/* 'readParameters' reads the parameters as entered in the command line. */
void readParameters(int argc, char* argv[]) {
	//for (int argumentIndex=1; argumentIndex<argc; argumentIndex++ ) { cout << argumentIndex  << ". " << argv[argumentIndex] << endl; }

	for (int argumentIndex=1; argumentIndex<argc; argumentIndex++ ) {
		string currentArgument = argv[ argumentIndex ];

		//find argument in parameterlist
		int parameterIndex = findParameter( currentArgument );
		if ( parameterIndex == -1 ) {
			cout << "unknown argument: " << currentArgument << endl;
			return;
		}

		if ( parameters[ parameterIndex ]->isUnary() ) {
			parameters[ parameterIndex ]->setValue(true); // default
			if ( (argumentIndex+1 < argc) && ( argv[ argumentIndex+1][0]!='-' ) ) { // so there are more arguments, and next one isn't regular -x
				if ( tolower(argv[ argumentIndex+1][0]) == 'f' || ( argv[ argumentIndex+1][0] =='0' ) ) {
					parameters[ parameterIndex ]->setValue(false);
				}
				argumentIndex++; // in any case increase the argument index
			}
		}
		else {   // argument needs a parameter
			argumentIndex++; // move on to next argument in the list
			if (argumentIndex >= argc ) {
				cout << "argument of " << currentArgument << " lacking.\n";
				return ;
			}
			if (argv[ argumentIndex ][0]=='-' ) {
				cout << "argument of " << currentArgument << " seems erroneous.\n";
				return ;
			}
			// but if everything is allright,
			//cout << "Giving " << currentArgument << " the value " << argv[ argumentIndex ] << endl;
			parameters[ parameterIndex ]->setValue( string(argv[ argumentIndex ]) );
		}
	}
}

/* 'printHelp' prints all parameters available. */
void printHelp() {
	cout << "\nProgram:    " << g_programName <<" (conversion of Pindel output to VCF format)\n";
	cout << "Version:    " << g_versionString << endl;
	cout << "Contact:    Eric-Wubbo Lameijer <e.m.w.lameijer@gmail.com>\n";
	cout << "Usage:      " << g_programName << " -p <pindel_output_file> -r <reference_file>\n";
	cout << "              -R <name_and_version_of_reference_genome> -d <date_of_reference_genome_version>\n";
	cout << "              [-v <vcf_output_file>]\n\n";
	cout << "           the -v parameter is optional; when no output file name is given, output is written\n";
	cout << "           to a file with the name <pindel_output_file>.vcf.\n\n";
	cout << "Example:    " << g_programName << " -p sample3chr20_D -r human_g1k_v36.fasta -R 1000GenomesPilot-NCBI36\n";
	cout << "              -d 20101123 -v sample3chr20_D.vcf\n\n";
	cout << "or (with -P): " << g_programName << " -P sample3chr20 -r human_g1k_v36.fasta -R 1000GenomesPilot-NCBI36\n";
	cout << "              -d 20101123 -v sample3chr20_all.vcf\n\n";
	cout << "Note:      -is only guaranteed to work correctly on output files produced by pindel version 0.2.3 and above.\n";
	cout << "           -LI and BP files (long insertion and break point files) have a different type of header and\n";
	cout << "            are not supported yet.\n\n";

	for (int parameterIndex=0; parameterIndex<parameters.size(); parameterIndex++ ) {
		parameters[ parameterIndex ]->describe();
	}
	exit( EXIT_SUCCESS );
}

bool isRegularPindelInput() {
	return parameters[ findParameter("-p") ]->isSet();
}
bool isRootPindelInput() {
	return parameters[ findParameter("-P") ]->isSet();
}

/* 'checkParameters' checks whether all required parameters have been set. */
bool checkParameters() {
	if (parameters[ findParameter("-h") ]->getBValue() ) {
		printHelp();
	}

	bool canRun = true;
	for (int parameterIndex=0; parameterIndex<parameters.size(); parameterIndex++ ) {
		if (parameters[ parameterIndex ]->isRequired() && !parameters[ parameterIndex ]->isSet()) {
			cout << "\nRequired parameter " << parameters[ parameterIndex ]->getShortName() << "/" << parameters[ parameterIndex ]->getLongName()
			     << " " << parameters[ parameterIndex ]->getDescription() << " needs to be set.\n\n";
			canRun = false;
		}  //if
	}

	if ( isRegularPindelInput() && isRootPindelInput() ) {
		cout << "Sorry, you can't use -p and -P at the same time, please choose one option.\n\n";
		canRun = false;
	}
	else if ( !isRegularPindelInput() && !isRootPindelInput() ) {
		cout << "Pindel2vcf needs a pindel input file, either use the -p or the -P option, please.\n\n";
		canRun = false;
	}

	if (!canRun) {
		cout << "For further information, please run " << g_programName <<" without arguments or with option -h/--help.\n\n";
	}
	return canRun;
}

/* 'setParameters' sets the filters to be used in the rest of the program. */
void setParameters() {
	g_par.vcffile = parameters[ findParameter( "-v" )]->getSValue();
	if (g_par.vcffile.compare("")==0) {
		if (isRegularPindelInput()) {
			g_par.vcffile = g_par.pindelfile + ".vcf";   // default
		}
		else if (isRootPindelInput()) {
			g_par.vcffile = g_par.pindelroot + ".vcf";
		}
		else {
			cout << "Error trying to construct output filename!\n";
			exit( EXIT_FAILURE );
		}
	}
}

/* 'throughFilter' checks whether the event is good enough to be written to the output file. */
bool throughFilter(SVData sv) {
	if (( g_par.minsize > 1 ) && ( abs( sv.getSize()) < g_par.minsize ) ) {
		return false;
	}
	if (( g_par.maxsize > 0 ) && ( abs( sv.getSize()) > g_par.maxsize ) ) {
		return false;
	}
	if ( g_par.bothstrands && !sv.bothStrands() ) {
		return false;
	}
	if ( ( g_par.minsuppSamples >= 1 ) && ( sv.getNumSupportSamples(g_par.onlyBalancedSamples, g_par.minimumStrandSupport) < g_par.minsuppSamples ) ) {
		return false;
	}
	if ( ( g_par.minsuppReads >= 1 ) && ( sv.getNumSupportReads() < g_par.minsuppReads ) ) {
		return false;
	}
	if ( ( g_par.maxSuppReads >= 1 ) && ( sv.getNumSupportReads() > g_par.maxSuppReads ) ) {
		return false;
	}
	if ( ( g_par.regionStart > 0 ) && ( sv.getPosition() < g_par.regionStart ) ) {
		return false;
	}
	if ( ( g_par.regionEnd > 0 ) && ( sv.getPosition() > g_par.regionEnd ) ) {
		return false;
	}
	if ( g_par.maxInterRepeatNo >= 0 && !sv.withinAllowedRepeatsInternal(g_par.maxInterRepeatLength, g_par.maxInterRepeatNo )) {
		return false;
	}
	if ( g_par.maxPostRepeatNo >= 0 && !sv.withinAllowedRepeatsPostIndel(g_par.maxPostRepeatLength, g_par.maxPostRepeatNo )) {
		return false;
	}

	// all filters passed
	return true;
}

/* 'makeSampleMap' converts a set of sample names to a map containing both sample names and the genotype of a sample. */
void makeSampleMap( const set<string>& sampleNames, map<string, int>& sampleMap ) {
	int count=0;
	for (set<string>::iterator setIt=sampleNames.begin(); setIt!=sampleNames.end(); setIt++ ) {
		sampleMap.insert( pair<string,int>( *setIt, count++) );
	}
}

void reportSVsInChromosome(
  const string& chromosomeID,
  const set<string>& chromosomeNames,
  const set<string>& sampleNames,
  InputReader& pindelInput,
  map< string, int >& sampleMap,
  Genome& genome,
  ofstream& vcfFile
) {
	// if no reads have been found for this chromosome, skip it
	if (chromosomeNames.find(chromosomeID) == chromosomeNames.end() ) {
		cout << "No reads for chromosome " << chromosomeID << ", skipping it.\n";
		return;
	}
	cout << "Processing chromosome " << chromosomeID << endl;
	// rewind file to start
	int regionStart = 0;
	int regionEnd = 0;
	SVData backupSV(sampleNames.size() );
	bool backupAvailable = false;
	do {
		cout << "reportSVsInChromosome: start reading region.\n";
		regionEnd = regionStart + g_par.windowSize*1000000;
		cout << "Reading region " << regionStart << "-" << regionEnd << endl;
		pindelInput.rewind();
		int counter=0;
		vector<SVData> svs;
		if (backupAvailable) {
			svs.push_back( backupSV );
		}
		while (!pindelInput.eof()) {
			SVData svd( sampleNames.size() ) ;
			convertIndelToSVdata( pindelInput, sampleMap, genome, svd, chromosomeID);
			if (!pindelInput.eof() && ( chromosomeID=="" || (svd.getChromosome()==chromosomeID && svd.getPosition()>=regionStart && svd.getPosition()<regionEnd)) ) {
				svs.push_back( svd );
			}
			counter++;
			//if (counter%10==0) cout << "At counter " << counter << " pos " << svd.getPosition() << endl;
		}

		cout << "Total reads: " << svs.size() << endl;
		sort ( svs.begin(), svs.end() );
		cout << "Sorting completed" << endl;
		// now output the SVs
		for (int svIndex=0; svIndex<svs.size(); svIndex++ ) {
			if ( svIndex!=svs.size()-1 && throughFilter( svs[ svIndex ]) ) {
				vcfFile << svs[ svIndex ];
			}
			else  {   //empty

			} // if else: whether the SV passes through the filter
		}  // for: loop over all SVs
		if (svs.size()>0) {
			backupSV = svs[ svs.size()-1 ];
			backupAvailable = true;
		}
		regionStart += (g_par.windowSize*1000000);
	}
	while (regionEnd<genome.getChromosome( chromosomeID )->size());
	if ( backupAvailable && throughFilter( backupSV) ) {
		vcfFile << backupSV;
	}
	//cout << "DEBUG:reportSVsInChromosome: exit.\n";
}

int main(int argc, char* argv[]) {
	initBaseArray();
	createParameters();
	readParameters(argc,argv);
	if (argc==1) {
		printHelp();
		exit(EXIT_SUCCESS);
	}
	if (!checkParameters()) {
		exit(EXIT_FAILURE);
	}
	setParameters();
	ofstream vcfFile(g_par.vcffile.c_str());
	set<string> sampleNames;
	set<string> chromosomeNames;

	InputReader pindelInput;
	if (isRegularPindelInput()) {
		pindelInput.addFile( g_par.pindelfile );
	}
	else if (isRootPindelInput()) {
		string rootFilename = g_par.pindelroot;
		pindelInput.addFile( rootFilename + "_D");
		pindelInput.addFile( rootFilename + "_SI");
		pindelInput.addFile( rootFilename + "_LI");
		pindelInput.addFile( rootFilename + "_INV");
		pindelInput.addFile( rootFilename + "_TD");
	}
	if (pindelInput.eof()) {
		cout << "The pindel file (-p) does not exist.\n";
		exit( EXIT_FAILURE );
	}
	//cout<< "Samples0:\n";
	//showSet( sampleNames );
	getSampleNamesAndChromosomeNames(pindelInput,sampleNames,chromosomeNames);
	cout<< "Samples:\n";
	showSet( sampleNames );
	cout << "Chromosomes in which SVs have been found:\n";
	showSet( chromosomeNames );

	map< string, int > sampleMap;
	makeSampleMap( sampleNames, sampleMap );
	createHeader(vcfFile,"pindel",g_par.referenceName, sampleNames);

	// read in reference
	Genome genome;
	readReference(g_par.reference,genome);

	if (g_par.chromosome != "" ) {
		// a specific chromosome has been specified
		reportSVsInChromosome( g_par.chromosome, chromosomeNames, sampleNames, pindelInput, sampleMap, genome, vcfFile );
	}

	else for (int chromosomeCount=0; chromosomeCount<genome.d_chromosomes.size(); chromosomeCount++ ) {
			reportSVsInChromosome( genome.d_chromosomes[ chromosomeCount ].getID(), chromosomeNames, sampleNames, pindelInput, sampleMap, genome, vcfFile );
			genome.d_chromosomes[ chromosomeCount ].removeFromMemory(); // to prevent memory overload
		}

	if (g_par.gatkCompatible) {
		cout << "\nNote: for this conversion, the -G (GATK) compatibility option was used; it's possible that this format is not compatible with other VCF-reading software. " <<
		     "Also note that GATK requires genotypes to be 0/0, 0/1 or 1/1 instead of undefined, like ./1 or . ('not detected'). However, since pindel cannot yet " <<
		     "genotype events " <<
		     "(distinguish between 0/1 and 1/1) all events are called as 0/0 (not found) or 0/1, even while some may very well be homozygous alternative (1/1).\n\n";
	}
	else {
		cout << "\nNote: for this conversion, the -G (GATK) compatibility option was not used; while this allows Pindel to indicate the uncertainty in genotypes, and should be " <<
		     "compatible with most software, this format " <<
		     "will not be compatible with GATK pipelines and tools such as GATK ValidateVariants; if you wish to input the vcf-file into the GATK pipeline, " <<
		     "please use the -G option.\n\n";
	}
}

