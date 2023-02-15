/*
 *  CChromosome.h
 *
 *  Created on: 2011-01-28
 *  Author: Alexandre Dionne-Laporte
 */

#ifndef CCHROMOSOME_H_
#define CCHROMOSOME_H_

class CVariant;
#include "CRunSettings.h"

#include <map>
#include <vector>
#include <string>
#include <iosfwd>
#include <iostream>
#include <sstream>
//#include <boost/thread/mutex.hpp>

using namespace std;

class EChromosome : public std::exception {
private:
	string _msg;
	string _file;
	string _line;
public:
	EChromosome(const string msg, const string file, const int line) : _msg(msg), _file(file) {
		std::stringstream ss;
		ss << line;
		_line = ss.str();
	};
	virtual ~EChromosome() throw() {}
//	virtual const char* what() const throw() {
//		string ret;
//		ret.append(_msg).append("\nIn file ").append(_file).append(", line ").append(_line);
//		return ret.c_str();
//	}
	const string message() {
		std::stringstream ss;
		ss << _msg << endl << "In file " << _file << ", line " << _line << endl;;
		return ss.str();
	}
};

class CChromosome {
	friend std::ostream & operator << (std::ostream &, const CChromosome &);

protected:
	string _name;
	multimap<int, CVariant*> * _variants;

private:
	static map<string, CChromosome*> * _chromosomes; /** for use by the singleton */
//	boost::mutex mtx_;

	/**
	 * Default constructor is private. To use with singleton
	 */
	CChromosome(string name) : _name(name) { _variants = new multimap<int, CVariant*>(); };
	~CChromosome() {
		_variants->clear();
		delete _variants;
//		if (CRunSettings::getInstance()->get_verbose_level() > CRunSettings::VERBOSE_LOG)
//			cerr << "CChromosome::~CChromosome" << endl;
	}

public:
	/** for use by the singleton */
	static CChromosome* getSetChromosome(string);
	static vector<CChromosome*> * getAllChromosomes() ;

	//GETTERS / SETTERS
	CVariant* getVariant(int position, string ref_allele, string novel_allele);
	CVariant* addVariant(CVariant*);
	CVariant* remove(CVariant*);

	const vector<CVariant*> * getAllVariants();
	inline string getName () const { return _name; };

};

#endif /* CCHROMOSOME_H_ */
