#include <string>
#include <ostream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath> 

#include "p2v.h"
#include "p2v_svdata.h"
#include "p2v_genome.h"

using std::string;
using std::endl;
using std::cout;
using std::min;
using std::max;  // genotype code
using std::sqrt; // genotype code
using std::stringstream; //genotype code

Genotype::Genotype() {
	Genotype( 0, 0, 0 );
}

Genotype::Genotype( const int readDepthPlus, const int readDepthMinus, const int totalRefSupport ) {
	d_readDepthPlus   = readDepthPlus;
	d_readDepthMinus  = readDepthMinus;
	d_totalRefSupport = totalRefSupport;
}

void Genotype::fuse( const Genotype& gt ) {
	d_readDepthPlus  += gt.d_readDepthPlus;
	d_readDepthMinus += gt.d_readDepthMinus;
	d_totalRefSupport = max(d_totalRefSupport, gt.d_totalRefSupport);
}

void Genotype::reset() {
		d_readDepthPlus   = 0;
		d_readDepthMinus  = 0;
		d_totalRefSupport = 0;
}

/** 'balanced' returns whether the a and b seem derived from a normal genomic binomial distribution of a and b, chance should be over 95% that the distribution is believable given heterozygosy (theoretically, as many
 reference strands as event strands. */
bool balanced( const unsigned int a, const unsigned int b) {
	if (a==b) {
		return true;
	}
	unsigned int smallest = min( a, b );
	unsigned int largest = max( a, b);
	unsigned int sum = a + b;
	if (sum>22) {
		return (smallest >= (sum/2) - sqrt( sum ));
	}
	else if (sum>19) {
		return ( smallest >= 6 );
	}
	else if (sum>16) {
		return ( smallest >= 5 );
	}
	else if (sum>14) {
		return ( smallest >= 4 );
	}
	else if (sum>11) {
		return ( smallest >= 3 );
	}
	else if (sum>8) {
		return ( smallest >= 2 );
	}
	else if (sum>5) {
		return ( smallest >= 1 );
	}
	else {
		return true;
	}
}

string deriveGenotype( const Genotype& rawGenotype )
// we may want to replace this simple logic by a function that assesses if reference support and event support are balanced, like 'isBalanced( eventSupp, refSupp );
{
	int totalEventReads = rawGenotype.getTotalReads();
	int totalReferenceReads = rawGenotype.getAverageRefSupport();

	if (totalEventReads + totalReferenceReads < getMinCoverage()) {
		return "0/0";
	}
	float AF = (float)totalEventReads / (totalEventReads + totalReferenceReads);
	if (AF < getHetCutoff()) {
		return "0/0";
	}
	else if (AF >= getHetCutoff() && AF < getHomCutoff()) {
		return "0/1";
	}
	else if (AF >= getHomCutoff()) {
		return "1/1";
	}

}

const string Genotype::getGTold() const {
	if (getGATKCompatible()) {
		if (d_readDepthPlus==0 && d_readDepthMinus==0) {
			return "0/0";
		}
		else {   // at least one read depth is 1
			return "0/1";
		}
	}
	else {
		if (d_readDepthPlus==0 && d_readDepthMinus==0) {
			return ".";
		}
		else {   // at least one genotype is 1
			return "1/.";
		}
	} // gatk compatibility not required

	return "ERROR"; // this should not happen
}

const string Genotype::getGTnew() const {
	return deriveGenotype( *this );
}

const string Genotype::getGTRDAD() const {
	stringstream ss;
	ss << getGTnew() << ":" << getTotalRefSupport() << "," << getTotalReads();
	return ss.str();
}

const string Genotype::getGTAD() const {
	stringstream ss;
	ss << getGTold() << ":" << getTotalReads();
	return ss.str();
}

string SVData::getAlternative() const {
	if (d_svtype == "INS" && d_svlen == 0)  { // long insertion
		return "<INS>";
	}
	string altVariant = "";
	const string* reference = d_genome_ptr->getChromosome( d_chromosome );
	if ( d_svtype == "INS" || d_svtype == "DEL" || d_svtype == "RPL" ) {
		if (!(getGATKCompatible() && altSameLengthAsRef())) {
			altVariant += (*reference)[ d_position ];
		}
		altVariant += d_nt;
	}
	else if ( d_svtype == "DUP:TANDEM" ) {
		string refVariant = getReference();
		altVariant += refVariant;
		altVariant += d_nt;
		altVariant += refVariant.substr(1); // skip the position before...
	}
	else if ( d_svtype == "INV" ) {
		string refVariant = getReference();
		if (getGATKCompatible() && altSameLengthAsRef()) {
			string complement = "";
			createComplement( refVariant, complement );
			altVariant = complement;
		}
		else {
			altVariant += (*reference)[ d_position ];
			altVariant += d_nt;
			string referenceInv = refVariant.substr(1);
			string complement = "";
			createComplement( referenceInv, complement );
			altVariant += complement;
			altVariant += d_nt2;
		}
	}
	return altVariant;
}


// getReference returns the reference sequence, for indels including the base before it
string SVData::getReference() const {
	const string* reference = d_genome_ptr->getChromosome( d_chromosome );
	if (d_svtype == "INS" && d_svlen == 0)  { // long insertion
		string refVariant = "";
		refVariant+= (*reference)[ d_position ];
		return refVariant;
	}
	else {   // normal insertion/deletion or whatever
		string refVariant="";
		int startPosition = d_position;
		if (getGATKCompatible() && altSameLengthAsRef() ) {
			startPosition = d_position + 1; // workaround GATK
		}
		for (int position=startPosition; position<d_end; position++ ) {
			refVariant += (*reference)[ position ];
		}
		return refVariant;
	}
}

void reportProblematicSVSize( int size, const string& chromosome, int position ) {
	if ( size >= getsizeToWarnFor() ) {
		cout << "Warning! The SV at chromosome " << chromosome << ", position " << position << " is of size " << size << ".";
		cout << " It won't be put in the VCF in base sequence detail, but in the format";
		cout << " 'chrom pos firstrefbase <SVType>'.";
		cout << " This behaviour can be overridden by setting the -co parameter to -1 or to a larger" ;
		cout << " value than the current SV size.\n";
	}
}

string SVData::getOutputFormattedReference() const {
	string defaultRef = getReference();
	string defaultAlt = getAlternative();
	reportProblematicSVSize( defaultRef.size(), getChromosome(), getPosition() );
	if (defaultAlt == "<INS>" ) {
		return defaultRef;
	}
	else {
		if (getCompactOutput()>1) {
			if ( defaultRef.size() > getCompactOutput() || defaultAlt.size() > getCompactOutput() ) {
				defaultRef.erase(1);
			}
		}
	}
	return defaultRef;
}

string SVData::getOutputFormattedAlternative() const {
	string defaultRef = getReference();
	string defaultAlt = getAlternative();

	reportProblematicSVSize( defaultAlt.size(), getChromosome(), getPosition() );
	if (defaultAlt == "<INS>" ) {
		return defaultAlt;
	}
	else {
		if ( getCompactOutput() > 1 ) {
			if ( defaultRef.size() > getCompactOutput() || defaultAlt.size() > getCompactOutput() ) {
				defaultAlt = "<" + d_svtype + ">";
			}
		}
	}
	return defaultAlt;
}

SVData::SVData(const int genotypeTotal) { // default settings
	d_id=".";
	d_quality=".";
	d_filter=".";
	d_replaceLen=0;
	d_homlen=0;
	d_end = 0;
	d_homseq="";
	d_nt="";
	d_nt2="";
	int numberOfSamples = genotypeTotal;
	if (numberOfSamples<=0) {
		numberOfSamples=1;
	}
	d_format.resize( numberOfSamples, Genotype(0,0,0) );
};


/* 'bothStrands' Is a SV supported by reads on both strands? */
bool SVData::bothStrands() const {
	bool strandPlus = false;
	bool strandMinus = false;

	for (int sampleIndex=0; sampleIndex<d_format.size(); sampleIndex++ ) {
		if ( d_format[ sampleIndex ].getReadDepthPlus() > 0 ) {
			strandPlus = true;
		}
		if ( d_format[ sampleIndex ].getReadDepthMinus() > 0 ) {
			strandMinus = true;
		}
	}
	return ( strandPlus && strandMinus );
}


/* 'getNumSupportSamples': how many samples support this SV? */
int SVData::getNumSupportSamples(const bool onlyBalancedSamples, const int minimumStrandSupport) const {
	int numSupportingSamples = 0;
	for (int sampleIndex=0; sampleIndex<d_format.size(); sampleIndex++ ) {
		int plusSupport = d_format[ sampleIndex ].getReadDepthPlus();
		int minSupport = d_format[ sampleIndex ].getReadDepthMinus();
		if (onlyBalancedSamples) {
			if (plusSupport>=minimumStrandSupport && minSupport>=minimumStrandSupport ) {
				numSupportingSamples++;
			}
		}
		else {   // Sample does not need to be balanced
			if (plusSupport>=minimumStrandSupport || minSupport>=minimumStrandSupport ) {
				numSupportingSamples++;
			}
		}
	}
	return numSupportingSamples;
}


/* 'getNumSupportReads': how many reads support this SV? */
int SVData::getNumSupportReads() const {
	int numReads = 0;
	for (int sampleIndex=0; sampleIndex<d_format.size(); sampleIndex++ ) {
		numReads += d_format[ sampleIndex ].getTotalReads();
	}
	return numReads;
}


/* comparison operator of SVData objects; helps sorting them for output. */
bool SVData::operator<(const SVData& otherSV ) const {
	if ( d_chromosome.compare( otherSV.d_chromosome ) != 0 ) {
		return ( d_chromosome.compare( otherSV.d_chromosome ) < 0 );
	}
	if ( getPosition() != otherSV.getPosition()) {
		return ( getPosition() < otherSV.getPosition() );
	}
	else {
		return ( d_svlen < otherSV.d_svlen );   // position equal: then do smallest SV first
	}
}

/* 'testHypothesis' tests whether the DNA string consists of a number of repeating "hypothesis" strings. If so, it
 returns the number of repeats; otherwise, it returns zero. */
int testHypothesis(const string& hypothesis, const string& sequence ) {
	int hypLen = hypothesis.size();
	int sequenceLen = sequence.size();

	for (int testedBaseIndex=0; testedBaseIndex<sequence.size(); testedBaseIndex++ ) {
		char currentHypothesisBase = hypothesis[ testedBaseIndex % hypLen ];
		if ( currentHypothesisBase != sequence[ testedBaseIndex ] ) {
			return 0;
		}
	}
	return sequenceLen/hypLen;
}


/* 'countRepeats' counts the repeats in a sequence. Returns the shortest result which can explain (allowing for rotation, CACAC => CA) the input sequence */
int countRepeats( const string& sequence, const int maxRepeatLength, int &bestSize ) {
	int maximumLen = min( maxRepeatLength, (int)(sequence.size()/2) );
	if (maxRepeatLength<0 ) {
		// no maximum repeat length set; assume it is infinite
		maximumLen = (int)(sequence.size()/2);
	}
	string hypothesis = "";
	int bestRepeatLength = 0;
	int bestRepeatNumber = 0;
	for (int repeatLen=1; repeatLen<=maximumLen; repeatLen++ ) {
		hypothesis += sequence[ repeatLen - 1 ];
		int repeats = testHypothesis( hypothesis, sequence );
		if (repeats>0 && repeats*hypothesis.size() > bestRepeatLength * bestRepeatNumber) {
			bestRepeatLength = hypothesis.size();
			bestRepeatNumber = repeats;
		}
	}
	bestSize = bestRepeatLength;
	return bestRepeatNumber;
}


/* 'getSVSequence' returns the inserted or deleted sequence; if the mutation is complex or weird, it returns the new ('alt') sequence. */
string SVData::getSVSequence() const {
	string ref = getReference();
	string alt = getAlternative();
	string modifiedSequence = "";
	int pos=0;
	int maxPos=min( ref.size(), alt.size() );
	while (pos < maxPos && ref[pos]==alt[pos]) {
		pos++;
	};
	if (pos == maxPos ) { // simple insertion or deletion...
		modifiedSequence = ( maxPos==ref.size() ? alt.substr(pos) : ref.substr(pos));
	}
	else {   // replacement
		modifiedSequence = alt.substr( pos );
	}
	return modifiedSequence;
}


/* 'withinAllowedRepeatsPostIndel' Is the number of repeats after the detected SV within the maximum allowed number of repeats (maxNoRepeats) of the
 fundamental repetitive unit of the SV? */
bool SVData::withinAllowedRepeatsPostIndel(const int maxRepeatLen, const int maxNoRepeats) const {
	// 1. get the minimum repeating unit from the inserted sequence
	string sequenceToAnalyze = getSVSequence();
	int actualRepeatLength = 0;
	int repeatCount = countRepeats( sequenceToAnalyze, maxRepeatLen, actualRepeatLength );
	if (actualRepeatLength>0) {
		string hypothesis = sequenceToAnalyze.substr( 0 , actualRepeatLength );
		int extendedRepeatCount = testHypothesis( hypothesis, sequenceToAnalyze + d_homseq );
		return (extendedRepeatCount - repeatCount <= maxNoRepeats );
	}
	else {
		int bestSize = 0;
		int extendedRepeatCount = countRepeats( sequenceToAnalyze + d_homseq, maxRepeatLen, bestSize );
		int repetitiveLength = bestSize * extendedRepeatCount;
		int repetitivePartPostIndel = repetitiveLength - sequenceToAnalyze.size();
		double noRepPostIndel = (double)repetitivePartPostIndel / bestSize;
		return (int)noRepPostIndel <= maxNoRepeats;
	}
}


/* 'withinAllowedRepeatsInternal' Is the number of repeats of the smallest repetitive unit in the detected SV within the maximum allowed number of repeats (maxNoRepeats)? */
bool SVData::withinAllowedRepeatsInternal(const int maxRepeatLen, const int maxNoRepeats) const {
	string sequenceToAnalyze = getSVSequence();
	int actualRepeatLength = 0;
	int repeatCount = countRepeats( sequenceToAnalyze, maxRepeatLen, actualRepeatLength );
	if ( repeatCount > maxNoRepeats ) {
		return false;
	}
	else {
		return true;
	}
}



/* first version of operator==. Are two events the same? We could make this more complicated, but first see if it works. */
bool SVData::operator==(const SVData& otherSV ) const {
	// for non-NT
	if ( ( d_svtype == "DEL" )
	     && ( otherSV.d_svtype == "DEL" )
	     && ( d_bpr_start == otherSV.d_bpr_start )
	     && ( d_bpr_end == otherSV.d_bpr_end )
	     && ( d_svlen == otherSV.d_svlen )
	     && ( d_chromosome.compare( otherSV.d_chromosome ) == 0 )
	     //( d_position == otherSV.d_position )
	     //&& ( d_chromosome.compare( otherSV.d_chromosome ) == 0 )
	     //&& ( d_svtype.compare( otherSV.d_svtype ) == 0 )
	     //&& ( d_end == otherSV.d_end )
	     //&& ( d_replaceLen == otherSV.d_replaceLen )
	   ) {
		return true;
	}
	// for NT
	if ( ( d_svtype == "RPL" )
	     && ( otherSV.d_svtype == "RPL" )
	     && ( ( d_svlen - d_replaceLen) == ( otherSV.d_svlen - otherSV.d_replaceLen ) )
	     && ( d_bpr_start == otherSV.d_bpr_start )
	     && ( d_chromosome.compare( otherSV.d_chromosome ) == 0 )
	   ) {
		return true;
	}
	if ( ( d_svtype == "INS" )
	     && ( otherSV.d_svtype == "INS" )
	     //&& ( ( d_svlen - d_replaceLen) == ( otherSV.d_svlen - otherSV.d_replaceLen ) )
	     && ( d_bpr_start == otherSV.d_bpr_start )
	     && ( d_bpr_end == otherSV.d_bpr_end )
	     && ( d_svlen == otherSV.d_svlen )
	     && ( d_chromosome.compare( otherSV.d_chromosome ) == 0 )
	   ) {
		return true;
	}
	return false;
}

/* 'fuse' takes over all parameters from the other (earlier-occurring and hence better) SV, but adds the total support. */
void SVData::fuse( SVData& otherSV ) {
	for (int sampleIndex=0; sampleIndex<d_format.size(); sampleIndex++ ) {
		otherSV.d_format[ sampleIndex ].fuse( d_format[ sampleIndex ] );
	}
	*this = otherSV;
}

int FACT(int n) {

	if (n == 0 || n == 1) {
		return 1;
	}
	int fact=1;
	int i;

	for (i=1; i<=n; i++) {
		fact*=i;
	}
	std::cout << n << " " << fact << std::endl;
	return fact;
}

double fisher_test(int a, int c, int b, int d) {
	double p = 0.0;
	int n = a + b + c + d;
	p = (FACT(a+b)*FACT(c+d)*FACT(a+c)*FACT(d+b)) / (double)(FACT(a)*FACT(b)*FACT(c)*FACT(d)*FACT(n));
	std::cout << a << " " << c << " " << b << " " << d << " " << p << std::endl;
	return p;
}


ostream& operator<<(ostream& os, const SVData& svd) {
	double somatic_p_value = 0.0;

	os << svd.d_chromosome << "\t";
	os << svd.getPosition() << "\t";
	os << svd.d_id << "\t";
	os << svd.getOutputFormattedReference() << "\t";
	os << svd.getOutputFormattedAlternative() << "\t";
	os << svd.d_quality << "\t";
	if (svd.d_format.size() == 2 && getSomatic()) {

		somatic_p_value = fisher_test(svd.d_format[0].getTotalReads(), svd.d_format[0].getTotalRefSupport(), svd.d_format[1].getTotalReads(), svd.d_format[1].getTotalRefSupport());
		//if (somatic_p_value < 0.05) svd.d_filter = "PASS";
	}
	if (somatic_p_value < 0.05) {
		os << "PASS\t";
	}
	else {
		os << svd.d_filter << "\t";
	}

	os << "END=" << svd.getVCFPrintEnd() << ";";
	os << "HOMLEN=" << svd.d_homlen << ";";
	if ( svd.d_homlen != 0 ) {
		os << "HOMSEQ=" << svd.d_homseq << ";";
	}

	os	<< "SVLEN=";
	if (svd.d_svtype.compare("RPL")==0 || svd.d_svtype.compare("DEL")==0 ) {
		if (svd.d_svlen > 0 ) {
			os << "-";  // deletions need negative SVLEN
		}
	}
	os << svd.d_svlen << ";";
	os << "SVTYPE=" << svd.d_svtype;
	if ( svd.d_svtype.compare("RPL")==0 || svd.d_svtype.compare("DUP:TANDEM")==0 || svd.d_svtype.compare("INV")==0 ) {
		os << ";NTLEN=" << svd.d_replaceLen;
	}
	if ( svd.d_svtype.compare("INV")==0 ) {
		os << "," << svd.d_replaceLenTwo;
	}

	if (svd.d_format.size() == 2 && getSomatic()) {
		os << ";" << somatic_p_value;
	}

	if (getpindel024uOrLater() && svd.getAlternative()!="<INS>") {
		os << "\tGT:AD";
	}
	else {
		os << "\tGT:AD";
	}

	for (int counter=0; counter<svd.d_format.size(); counter++ ) {
		os << "\t";
		if (getpindel024uOrLater() && svd.getAlternative()!="<INS>") {
			os << svd.d_format[ counter ].getGTRDAD();
		}
		else {
			os << svd.d_format[ counter ].getGTAD();
		}
	}

	os << endl;

	return os;
}



