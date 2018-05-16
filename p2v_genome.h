#ifndef P2V_GENOME_H
#define P2V_GENOME_H

#include <string>
#include <vector>

#include "p2v_chromosome.h"


using std::string;
using std::vector;

/* 'Genome' contains a collection of chromosomes, and returns a pointer to the requested chromosome (the chromosome with the requested ID) */
class Genome {
public:
	const string* getChromosome( const string& id );
	void addChromosome( const Chromosome& chr ) {
		d_chromosomes.push_back( chr );
	}
	string firstChromosomeName();
	string nextChromosomeName();

	vector<Chromosome> d_chromosomes; // working fast here, but ideally I'd keep this private and implement an iterator
};

#endif
