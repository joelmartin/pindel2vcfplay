#ifndef P2V_CHROMOSOME_H
#define P2V_CHROMOSOME_H

#include <string>
using std::string;

/* 'Chromosome' contains a string identifier as well as the base sequence of the chromosome itself. */
class Chromosome {
public:
	Chromosome(const string& identifier, const string& fastaFilename, const int& length, const int& offset) {
		d_identifier = identifier;
		d_sequence=NULL;
		d_fastaFilename=fastaFilename;
		d_offset = offset;
		d_length = length;
	}
	~Chromosome() {
		delete d_sequence;
	}
	const string* getChromPtr();
	const string getID() const {
		return d_identifier;
	}

	void removeFromMemory();

private:
	void readFromFile();
	void readFromFile_test();
	string d_identifier;
	string* d_sequence;
	string d_fastaFilename;
  int d_offset;
	int d_length;
};

void initBaseArray();
void makeStrangeBasesN(string& dna);
bool normalBase(char ch );

#endif
