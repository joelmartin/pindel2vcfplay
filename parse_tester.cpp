#include <sstream>
#include <vector>
#include <string> 
#include <iostream>
#include <fstream>
#include <map>
#include <set>

#include "GetPot"
#include "p2v.h"
#include "p2v_input_reader.h"

using std::string;
using std::ios_base;
using std::ifstream;
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

void read_fai( map<string, int> &fa_map, const string fai_name );
bool isPindelSVIdentifier( string identifier );
bool isSVSummarizingLine( string line );
string fetchElement( istream& instream, const int index );
int countElements( istream& instream );
void parse_to_just_type(InputReader& pindelInput, set<string>& sampleNames, set<string>& chromosomeNames );

const int FIRST_SAMPLE_INDEX = 32;

int PARSE_BEEN_CALLED = 0;


int main (int argc, char **argv) {
  ios_base::sync_with_stdio(false);

	GetPot   cl(argc, argv);

	// fasta=whatever.fa
	// ctg=scaffold_8366
    // pindel=/Users/j_martin/?
    
	ifstream index_file, fasta_file;
	map< string, int > fa_map;

	const string fa  = cl("fasta", "default");
	const string ctg_id = cl( "ctg", "" );
	const string fai = fa + ".fai";
    const int    jobs = cl("jobs",10000);
    
	read_fai( fa_map, fai );
	map<string,int>::iterator it = fa_map.find( ctg_id );
	if ( it == fa_map.end() ) {
		cerr << "not found: " << ctg_id << "\n";
		exit( 1 );
	}
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
    //cout<< "Samples0:\n";
    //showSet( sampleNames );
    for( int i = 0; i< jobs; i++)
        parse_to_just_type(pindelInput,sampleNames,chromosomeNames);
    cout << "done: " << PARSE_BEEN_CALLED << endl;

  return 0;
}

void parse_to_just_type(InputReader& pindelInput, set<string>& sampleNames, set<string>& chromosomeNames ) {
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
// read fai and fill map

void read_fai( map<string, int> &fa_map, const string fai_name ) {
	struct tokens {
		string ctg;
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
            ss >> buffy.ctg, buffy.offset, buffy.offset;
			fa_map.insert( make_pair( buffy.ctg, buffy.offset ) );
		}
	}
}

