#include <string>
using std::string;

/* 'Chromosome' contains a string identifier as well as the base sequence of the chromosome itself. */
class Chromosome {
public:
	Chromosome(const string& identifier, const string& fastaFilename) {
		d_identifier = identifier;
		d_sequence=NULL;
		d_fastaFilename=fastaFilename;
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
	string d_identifier;
	string* d_sequence;
	string d_fastaFilename;
};

void initBaseArray();
void makeStrangeBasesN(string& dna);
bool normalBase(char ch );

