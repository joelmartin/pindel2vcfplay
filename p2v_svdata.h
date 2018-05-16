#ifndef P2V_SVDATA_H
#define P2V_SVDATA_H

#include <string>
#include <sstream> 
#include <ostream>

#include "p2v.h"
//#include "p2v_genotype.h"
#include "p2v_genome.h"

using std::ostream;
using std::ostringstream;
using std::string;

/* 'Genotype' stores the genotype of a mutation and its read depth. So, for example, "0/1:5". */
class Genotype {

	friend ostream& operator<<(ostream& os, const Genotype& gt);

public:
	Genotype();
	Genotype( const int readDepthPlus, const int readDepthMinus, const int totalRefSupport );
	void fuse( const Genotype& gt );
	void reset();
	int getReadDepthPlus() const {
		return d_readDepthPlus;
	}
	int getReadDepthMinus() const {
		return d_readDepthMinus;
	}
	int getTotalReads() const {
		return ( d_readDepthPlus + d_readDepthMinus ) ;
	}

	int getAverageRefSupport() const {
		return d_totalRefSupport;
	}

	int getTotalRefSupport() const {
		return d_totalRefSupport;
	}
	const string getGTRDAD() const;
	const string getGTAD() const;

private:
	int d_readDepthPlus;
	int d_readDepthMinus;
	int d_totalRefSupport;
	const string getGTold() const;
	const string getGTnew() const;
};

/* 'SVData' stores the data of a certain structural variant. */
class SVData {

	friend ostream& operator<<(ostream& os, const SVData& svd);

public:
	SVData(const int genotypeTotal);

	int getPosition() const {
		// attempted workaround for GATK undocumented feature
		if (getGATKCompatible() && altSameLengthAsRef()) {
			return d_position+1;
		}
		else {
			return d_position;
		}
	}
	int getSize() const {
		return d_svlen;
	}
	bool bothStrands() const;
	int getNumSupportSamples(const bool onlyBalancedSamples, const int minimumStrandSupport) const;
	int getNumSupportReads() const;

	string getChromosome() const {
		return d_chromosome;
	}

	void setChromosome(const string& chromosome) {
		d_chromosome = chromosome;
	}
	void setPosition(const int position) {
		d_position = position;
	}
	void setID(const string& id) {
		d_id = id;
	}

	void setQuality(const double quality) {
		ostringstream strs;
		strs << quality;
		d_quality = strs.str();
	}
	void setFilter(const string& filter) {
		d_filter = filter;
	}

	void setEnd(const int end) {
		//cout << "Setting end to " << end << endl;

		d_end = end;
	}
	int getVCFPrintEnd() const {
		if (d_end <= d_position ) {
			//cout << "Warning: end position of the SV (" << d_end << ") appears to be before the startposition (" << d_position << "). Will adjust end to be startposition+reflen-1.\n";
		}
		else {}   // empty else.
		return int( d_position+getReference().size()-1 ); //explicitize.jam
	}

	void setHomlen(const int homlen) {
		d_homlen = homlen;
	}
	void setHomseq(const string& homseq) {
		d_homseq = homseq;
	}
	void setSVlen(const int svlen) {
		d_svlen = svlen;
	}
	void setSVlen(const string& svlen) {
		d_svlen = atoi(svlen.c_str());
	}
	void setSVtype(const string& svtype) {
		d_svtype = svtype;
	}

	//void setSupportingReads(const int supportingreads) { d_supportingreads = supportingreads; }
	void addGenotype(const int sampleID, const int readDepthPlus, const int readDepthMinus, const int totalRefSupport) {
		Genotype gt(readDepthPlus,readDepthMinus, totalRefSupport);
		d_format[ sampleID ] = gt ;
	}

	void setReplace( int replaceLength, int secondReplaceLen=-1 ) {
		d_replaceLen = replaceLength;
		d_replaceLenTwo = secondReplaceLen;
	}
	void setBPrange( const int bpr_start, const int bpr_end ) {
		d_bpr_start = bpr_start;
		d_bpr_end = bpr_end ;
	}
	void setGenome(Genome& genome) {
		d_genome_ptr = &genome;
	};
	bool operator<(const SVData& otherSV ) const;
	bool operator==(const SVData& otherSV ) const;
	void fuse( SVData& otherSV );
	string getReference() const; // reference sequence, for indels including the base before it
	string getAlternative() const;
	string getOutputFormattedReference() const;
	string getOutputFormattedAlternative() const;
	void setNT(string nt) {
		d_nt = nt;
	};
	void setSecondNT(string secondNT) {
		d_nt2 = secondNT;
	};
	bool withinAllowedRepeatsPostIndel(const int maxRepeatLen, const int maxNoRepeats) const;
	bool withinAllowedRepeatsInternal(const int maxRepeatLen, const int maxNoRepeats) const;
	bool isEquilengthReplacement() const {
		if ((d_svtype=="RPL" && d_svlen == d_replaceLen ) || (d_svtype=="INV" && d_replaceLen==0 && d_replaceLenTwo==0 )) {
			return true;
		}
		else {
			return false;
		}
	}

private:

	bool altSameLengthAsRef() const {
		return ((d_svtype=="RPL" && d_svlen == d_replaceLen) ||
		        (d_svtype=="INV" && d_replaceLen==0 && d_replaceLenTwo ==0 ));
	}

	int d_position; // 1-based

	int d_end;
	int d_homlen;
	int d_bpr_start, d_bpr_end;
	int d_svlen;
	int d_replaceLen;
	int d_replaceLenTwo; // for inversions
	double p_somatic;

	string d_nt;
	string d_nt2; // for inversions
	vector<Genotype> d_format;
	Genome* d_genome_ptr;
	string d_chromosome;
	string d_id; // default '.', as we don't mine variant databases yet
	string d_quality; // '.' by default, but can be floating-point number
	string d_filter;  // "PASS" by default
	string d_svtype;
	string d_homseq;

	string getSVSequence() const;
};

#endif
