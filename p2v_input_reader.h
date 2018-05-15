#ifndef INPUT_READER_H
#define INPUT_READER_H

#include <string>
#include <vector>
#include <fstream>

using std::string;
using std::vector;
using std::ifstream;

/** 'InputReader' can house a vector of files, allowing access as if it were one huge file. */
class InputReader {
public:
	InputReader();
	string getLine();
	bool eof();
	void addFile(const string filename);
	void rewind();

private:
	vector<string> m_filenames;
	int m_nextFileIndex;
	bool m_readable;
	ifstream m_currentFile;
	bool canReadMore();
	void moveToNextFile();
};

#endif
