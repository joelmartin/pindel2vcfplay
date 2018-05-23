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
	//profiling cruft
  int md_moved_to_next_file = 0;
  int md_rewound = 0;
  int md_opened = 0;
  int md_gotline = 0;
    
private:
	vector<string> m_filenames;
	vector<int>    m_positions;
	int m_nextFileIndex;
	bool m_readable;
	ifstream m_currentFile;
	bool canReadMore();
	void moveToNextFile();

};

#endif
