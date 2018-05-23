#include <string>
#include <iostream>
#include <fstream>
//#include <sstream>
#include <vector>

#include "p2v_input_reader.h"

using std::string;
using std::vector;
using std::ifstream;

string InputReader::getLine() {
    md_gotline++;
	if (canReadMore()) {
		string line;
		getline( m_currentFile, line );
		return line;
	}
	else {
		return "";
	}
}

bool InputReader::canReadMore() {
	// default case: current file is okay
	if (m_currentFile && !m_currentFile.eof()) {
		return true;
	}

	while (!m_currentFile || m_currentFile.eof()) {
		moveToNextFile();
		if (!m_readable) {
			break;   // final EOF
		}
	}

	return m_readable;
}

void InputReader::moveToNextFile() {
    md_moved_to_next_file++;
	if (m_nextFileIndex<m_filenames.size()) {
		m_currentFile.close();
		m_currentFile.open( m_filenames[ m_nextFileIndex ].c_str() );
		m_nextFileIndex++;
    md_opened++;
	}
	else {
		m_readable = false;
	}
}

void InputReader::rewind() {
	m_currentFile.open("");
	m_nextFileIndex = 0;
	m_readable = true;
	md_rewound++;
}

InputReader::InputReader() {
	rewind();
}

void InputReader::addFile(const string filename) {
	m_filenames.push_back( filename );
}

bool InputReader::eof() {
	return !canReadMore();
}

