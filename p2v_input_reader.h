#ifndef INPUT_READER_H
#define INPUT_READER_H

#include <string>
#include <vector>
#include <fstream>
#include <map>
//#include <pair>

using std::string;
using std::vector;
using std::map;
using std::pair;
using std::ifstream;

typedef std::map<string, int> chrPos;  // Contig_Name, offset
typedef std::vector<chrPos>  filemaps;

/** 'InputReader' can house a vector of files, allowing access as if it were one huge file. */
class InputReader {
public:
	InputReader();
	string getLine();
	bool eof();
	void addFile(const string filename);
	void addChrPos( string chromosomeName );
	void rewind();
	bool chrSeenInFile( string chromosomeName ); // is it in file.filemaps
	void addFileChrPos( string chromosomeName ); // not used? addChrPos
	void setChrTarget( string chromosomeName );
	void initFiles();
	void pastCID();
	void setCID(string chromosomeName);  //unused but maybe better name thatn setChrTarget
	//profiling cruft
	int md_moved_to_next_file = 0;
	int md_rewound = 0;
	int md_opened = 0;
	int md_gotline = 0;
	
private:
	vector<string> m_filenames;
	filemaps m_positions;
	int m_nextFileIndex;
	int m_currentFileIndex;
	int m_prevPos;               // position before current line was read.
	bool m_readable;
	ifstream m_currentFile;
	bool canReadMore();
	void moveToNextFile();
	string m_chromosomeTarget;
};

#endif
