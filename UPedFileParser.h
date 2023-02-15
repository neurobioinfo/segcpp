/*
 * CPedFile.h
 *
 *  Created on: 2011-10-31
 *  Author: Alexandre Dionne-Laporte
 */

#ifndef CPEDFILE_H_
#define CPEDFILE_H_

class CSample;

#include <fstream>
#include <exception>
#include <vector>

using namespace std;

class EPedFileParser : public std::exception {
private:
	string _msg;
public:
//	virtual const char* what() const throw() {
//		return _msg.c_str();
//	}
	const string message() {
		std::stringstream ss;
		ss << _msg ;
		return ss.str();
	}
	EPedFileParser(const string msg) : _msg(msg) {};
	virtual ~EPedFileParser() throw() {}
};

class UPedFileParser {
public:
	static vector<CSample*> * parseFile ( std::ifstream& inStream) throw (EPedFileParser) ;
};

#endif /* CPEDFILE_H_ */
