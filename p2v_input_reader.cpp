#include <string>
#include <iostream>
#include <fstream>
//#include <sstream>
#include <vector>
#include <algorithm>

#include "p2v_input_reader.h"

using std::string;
using std::vector;
using std::ifstream;
using std::make_pair;
//using std::pair;
using std::find;

string InputReader::getLine() {
	md_gotline++;
	if (canReadMore()) {
		string line;
		m_prevPos = m_currentFile.tellg();
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
void InputReader::pastCID(){  //having caller call this when past chromosome of interest.
	moveToNextFile();
}
void InputReader::moveToNextFile() {
	md_moved_to_next_file++;
	if (m_nextFileIndex<m_filenames.size()) {
		m_currentFile.close();
		m_currentFileIndex = m_nextFileIndex;
		m_currentFile.open( m_filenames[ m_nextFileIndex ].c_str() );
		md_opened++;
		m_nextFileIndex++;
		
		if ( ! m_chromosomeTarget.empty() ) {
			if ( ! m_positions.empty() && ! m_positions[ m_currentFileIndex ].empty() ) {
				auto it = m_positions[ m_currentFileIndex ].find( m_chromosomeTarget );
				if ( it != m_positions[ m_currentFileIndex].end() ) {
					m_currentFile.seekg( it->second );
				}
			}
		}
		else {
			m_prevPos = 0;
		}
	}
	else {
		m_readable = false;
	}
}

bool InputReader::chrSeenInFile( string chromosomeName ) {
	bool foundit = false;
	if ( ! m_positions.empty() && ! m_positions[ m_currentFileIndex ].empty() ) {
		auto it = m_positions[ m_currentFileIndex ].find( chromosomeName );
		if ( it != m_positions[ m_currentFileIndex].end() ) {
			foundit = true;
		}
		else {
			foundit = false;
		}
	}
	else {
		foundit =  false;
	}
	return( foundit );
}
void InputReader::setChrTarget( string chromosomeName ) {
	m_chromosomeTarget = chromosomeName;
}
void InputReader::rewind() {
	m_currentFile.open("");  // starts out open now
	m_nextFileIndex = 0;
	m_currentFileIndex = 0;
	m_prevPos = 0;
	m_readable = true;
	md_rewound++;
}

void InputReader::addChrPos( string chromosomeName) {
	m_positions[ m_currentFileIndex].insert( make_pair(chromosomeName, m_prevPos) );
}

InputReader::InputReader() {
	rewind();
}

void InputReader::addFile(const string filename) {
	chrPos tempPos;
	m_filenames.push_back( filename );
	//m_files.push_back(make_shared<ifstream>( filename ) );
	m_positions.push_back(tempPos);
}

bool InputReader::eof() {
	return !canReadMore();
}

