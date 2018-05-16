
#include <string>
#include "p2v_genome.h"

using std::string;

const string* Genome::getChromosome( const string& id ) {

	for (int chromosomeIndex=0; chromosomeIndex<d_chromosomes.size(); chromosomeIndex++ ) {
		if ( id.compare( d_chromosomes[ chromosomeIndex ].getID() ) == 0 ) {
			return d_chromosomes[ chromosomeIndex ].getChromPtr();
		}
	}
	return NULL; // default value if chromosome not found
}


