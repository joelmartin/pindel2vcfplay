#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <set>

#include "GetPot"
//#include "p2v.h"
#include "p2v_input_reader.h"

using std::string;
using std::ios_base;
using std::ifstream;
using std::ofstream;
using std::istream;
using std::ostream;
using std::cin;
using std::cout;
using std::cerr;
using std::map;
using std::set;
using std::stringstream;
using std::vector;
using std::make_pair;
using std::endl;

struct fai_s {
	int offset;
	int length;
};

void read_fai( map<string, int> &fa_map, const string fai_name );
void read_fai2( map<string, fai_s> &fa_map, vector <string> &allChromosomeNames, const string fai_name );
bool isPindelSVIdentifier( string identifier );
bool isSVSummarizingLine( string line );
string fetchElement( istream& instream, const int index );
int countElements( istream& instream );
void GetSampleNamesAndChromosomeNames(InputReader& pindelInput, set<string>& sampleNames, set<string>& chromosomeNames );
void convertIndelToSVdata( InputReader& pindelInput, const string& targetChromosomeID, ofstream& debug);
void reportSVsInChromosome( const string& chromosomeID, const set<string>& chromosomeNames, const set<string>& sampleNames,InputReader& pindelInput, map<string,fai_s>& fai_map, ofstream& debug);
//rootFilename, chromosomeNames, sampleNames, pindelInput, fai_map );
const int FIRST_SAMPLE_INDEX = 32;

int PARSE_BEEN_CALLED = 0;

/*struct fai {
 fai() : offset(0), length(0) {}
 fai(int newOffset, int newLength)
 : offset(newOffset), length(newLength) {}
 
 int offset;
 int length;
 };*/


//int main() {
//categories[1] = category(1, "First category");
//	categories[2] = category(2, "Second category");

int main (int argc, char **argv) {
	ios_base::sync_with_stdio(false);
	
	GetPot   cl(argc, argv);
	
	// fasta=whatever.fa
	// ctg=scaffold_8366
	// pindel=/Users/j_martin/?
	string debug_name = "debug.out";
	ofstream debug(debug_name);
	if ( ! debug.good() ) {
		cerr << "failed to open: " << debug_name;
		exit( 1 );
	}
	
	map< string, int > fa_map;
	map< string, fai_s > fai_map;
	vector <string> allChromosomeNames;
	
	const string fa  = cl("fasta", "default");
	const string ctg_id = cl( "ctg", "" );
	const string fai = fa + ".fai";
	const int    jobs = cl("jobs",10000);
	
	read_fai( fa_map, fai );
	read_fai2( fai_map, allChromosomeNames, fai);
	map<string,int>::iterator it = fa_map.find( ctg_id );
	map<string,fai_s>::iterator it2 = fai_map.find(ctg_id);
	if ( it2 == fai_map.end()){
		cerr << "not found: " << ctg_id << "\n";
		exit( 1 );
	}
	if ( it == fa_map.end() ) {
		cerr << "not found: " << ctg_id << "\n";
		exit( 1 );
	}
	cerr << "found: " << ctg_id << " length: " << it2->second.length << " offset: " << it2->second.offset << endl;
	/* testing parsing */
	const string rootFilename = cl("pindel","");
	InputReader pindelInput;
	set<string> sampleNames;
	set<string> chromosomeNames;
	
	pindelInput.addFile( rootFilename + "_D");
	pindelInput.addFile( rootFilename + "_SI");
	pindelInput.addFile( rootFilename + "_LI");
	pindelInput.addFile( rootFilename + "_INV");
	pindelInput.addFile( rootFilename + "_TD");
	if (pindelInput.eof()) {
		cout << "The pindel file (-p) does not exist.\n";
		exit( EXIT_FAILURE );
	}
	GetSampleNamesAndChromosomeNames(pindelInput,sampleNames,chromosomeNames); //parse pindelInput 1
	cerr << "parse called: " << PARSE_BEEN_CALLED << endl;
	cerr << "movedtos: " << pindelInput.md_moved_to_next_file << endl;
	cerr << "rewounds: " << pindelInput.md_rewound << endl;
	cerr << "openeds:  " << pindelInput.md_opened << endl;
	cerr << "gotlines: " << pindelInput.md_gotline << endl;
	
	/*	for ( set<string>::iterator it=chromosomeNames.begin(); it!=chromosomeNames.end(); ++it ) {  //oops ordered
	 reportSVsInChromosome( *it, chromosomeNames, sampleNames, pindelInput, fai_map, debug );
	 }
	 */
	for ( vector<string>::iterator it=allChromosomeNames.begin(); it!=allChromosomeNames.end(); ++it ) {  //oops ordered
		reportSVsInChromosome( *it, chromosomeNames, sampleNames, pindelInput, fai_map, debug );
	}
	cerr << "parse called: " << PARSE_BEEN_CALLED << endl;
	cerr << "movedtos: " << pindelInput.md_moved_to_next_file << endl;
	cerr << "rewounds: " << pindelInput.md_rewound << endl;
	cerr << "openeds:  " << pindelInput.md_opened << endl;
	cerr << "gotlines: " << pindelInput.md_gotline << endl;
	
	return 0;
}

void convertIndelToSVdata( InputReader& pindelInput, const string& targetChromosomeID, ofstream& debug) {
	string line;
	do {
		line = pindelInput.getLine();
	}
	while (!pindelInput.eof() && (! isdigit(line[0]) || !isSVSummarizingLine( line )));
	
	if (pindelInput.eof()) {
		return;
	}
	// at this point have a summarizing line
	// in a perfect world file will end up with one line for every item.
	debug << line << endl;
	// does a lot of stuff with that line here
	stringstream lineStream;
	lineStream << line;
	string svType = fetchElement( lineStream, 2 ); // to 2
	if ( svType.compare("LI") == 0 ) {
		string chromosomeID = fetchElement( lineStream, 2);
		if (chromosomeID!=targetChromosomeID) {
			pindelInput.pastCID();
			return;
		}
		else {
			debug << "evaluating " << line << endl;
		}
		int beforeStartPos = atoi( fetchElement( lineStream, 1 ).c_str() );
		cout << "hey, found LI line" << endl;
		return;
	}
	string svLen = fetchElement( lineStream, 1 ); //jam.is that int-able?
	string numNTaddedStr = fetchElement( lineStream, 2 );
	string ntAdded = fetchElement( lineStream, 1 ); // to 6
	string chromosomeID = fetchElement( lineStream, 2); // now at position 8
	if (chromosomeID!=targetChromosomeID) {
		pindelInput.pastCID();
		return;
	}
	else {
		debug << "evaluating " << line << endl;
	}
}


void reportSVsInChromosome( const string& chromosomeID,
													 const set<string>& chromosomeNames,
													 const set<string>& sampleNames,
													 InputReader& pindelInput,
													 map<string,fai_s>& fai_map,
													 ofstream& debug
													 ){
	// if no reads have been found for this chromosome, skip it
	if (chromosomeNames.find(chromosomeID) == chromosomeNames.end() ) {
		cout << "No reads for chromosome " << chromosomeID << ", skipping it.\n";
		return;
	}
	cout << "Processing chromosome " << chromosomeID << endl;
	// rewind file to start
	int regionStart = 0;
	int regionEnd = 0;
	map<string,fai_s>::iterator fai_map_it = fai_map.find( chromosomeID );
	pindelInput.setChrTarget( chromosomeID);
	do {
		cout << "reportSVsInChromosome: start reading region.\n";
		regionEnd = regionStart + 300*1000000;
		cout << "Reading region " << regionStart << "-" << regionEnd << endl;
		pindelInput.rewind(); // change to go to last position read?
		int counter=0;
		while (!pindelInput.eof()) {
			convertIndelToSVdata( pindelInput, chromosomeID, debug);
			counter++;
			//if (counter%10==0) cout << "At counter " << counter << " pos " << svd.getPosition() << endl;
		}
		regionStart += (300*1000000);
		//} while ( !regionEnd<genome.getChromosome( chromosomeID )->size());
	} while ( regionEnd < fai_map_it->second.length );
}

void GetSampleNamesAndChromosomeNames(InputReader& pindelInput, set<string>& sampleNames, set<string>& chromosomeNames ) {
	//cout << "DEBUG:start GetSampleNamesAndChromosomeNames\n";
	string line;
	unsigned int counter=0;
	
	// loop over the input file, registering each 'interesting' pindel output line (pindel output has the structure (ABCDn)m
	// A is a ##### separator line, B the interesting line (with sample names and coverage), C the line with the reference DNA
	// sequence (starting with ACTG) and the D-lines start with spaces and give the individual reads that support an event
	// HOWEVER, it is also possible that an user greps only the 'useful' output lines.
	
	PARSE_BEEN_CALLED++;
	while (!pindelInput.eof()) {
		
		// find the next 'B-line' (the line containing the sample names, coverage etc.)
		do {
			line = pindelInput.getLine();
			//cout << "isdigit: " << isdigit(line[0]) << " ,line: " << line;
		}
		while (!pindelInput.eof() && (! isdigit(line[0]) || !isSVSummarizingLine( line )));
		
		if (pindelInput.eof()) {
			//cout << "DEBUG:end GetSampleNamesAndChromosomeNames\n";
			return;
		}
		
		// 'line' should contain a B-line now. Analyze it.
		counter++;
		//cout << "C: " << counter << endl;
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
		if ( ! pindelInput.chrSeenInFile(chromosomeName) ) {
			pindelInput.addChrPos(chromosomeName);
		}
		chromosomeNames.insert( chromosomeName );
		//		cout << "Studying chromosome " << chromosome << endl;
		
		// 8 = 2+6, so corrects for previous reads
		
		
		int numberOfSamples = atoi( fetchElement( lineStream, FIRST_SAMPLE_INDEX - 12 ).c_str() );
		string firstSampleName = fetchElement( lineStream, 4 );
		if ( firstSampleName != "" ) {
			sampleNames.insert( firstSampleName );
		}
		
		//cout << "ElInLine: " << elementsInLine << ", FSINDEX: " << FIRST_SAMPLE_INDEX << ", NoS=" << numberOfSamples << endl;
		bool pindel024uOrLater = false;
		
		if (elementsInLine> FIRST_SAMPLE_INDEX + 5 * numberOfSamples ) {
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
	return;
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
void read_fai2( map<string, fai_s> &fa_map, vector<string>& allChromosomeNames, const string fai_name ) {
	struct tokens {
		string ctg;
		int length;
		int offset;
	};
	fai_s fai_struct;
	tokens buffy;
	ifstream index_file;
	string line;
	
	index_file.open( fai_name, ifstream::in );
	
	if ( ! index_file.good() ) {
		cerr << "failed to open: " << fai_name << "\n";
		exit(1);
	}
	else {
		while ( getline( index_file, line ) ) {
			stringstream ss(line);
			ss >> buffy.ctg;
			ss >> fai_struct.length;
			ss >> fai_struct.offset;
		
			allChromosomeNames.push_back(buffy.ctg);
			fa_map.insert( make_pair( buffy.ctg, fai_struct ) );
		}
	}
}

// read fai and fill map

void read_fai( map<string, int> &fa_map, const string fai_name ) {
	struct tokens {
		string ctg;
		int length;
		int offset;
	};
	tokens buffy;
	ifstream index_file;
	string line;
	
	index_file.open( fai_name, ifstream::in );
	
	if ( ! index_file.good() ) {
		cerr << "failed to open: " << fai_name << "\n";
		exit(1);
	}
	else {
		while ( getline( index_file, line ) ) {
			stringstream ss(line);
			ss >> buffy.ctg;
			ss >> buffy.length;
			ss >> buffy.offset;
			fa_map.insert( make_pair( buffy.ctg, buffy.offset ) );
		}
	}
}

