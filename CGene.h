/*
 *  CGene.h
 *
 *  Created on: 2011-01-28
 *  Author: Alexandre Dionne-Laporte
 */

#ifndef CGENE_H_
#define CGENE_H_

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

class EGene : public std::exception {
private:
	string _msg;
	string _file;
	string _line;
public:
	EGene(const string msg, const string file, const int line) : _msg(msg), _file(file) {
		std::stringstream ss;
		ss << line;
		_line = ss.str();
	};
	virtual ~EGene() throw() {}
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

class CGene {
	friend std::ostream & operator << (std::ostream &, const CGene &);

protected:
	string _name;
	vector<CVariant*> * _variants;

private:
	static map<string, CGene*> * _genes; /** for use by the singleton */
//	boost::mutex mtx_;

	/**
	 * Default constructor is private. To use with singleton
	 */
	CGene(string name) : _name(name) {
		_variants = new vector<CVariant*>();
	};

	~CGene() {
		_variants->clear();
		delete _variants;
//		if (CRunSettings::getInstance()->get_verbose_level() > CRunSettings::VERBOSE_LOG)
//			cerr << "CGene::~CGene" << endl;
	};

public:
	/** for use by the singleton */
	static CGene* getSetGene(string);
	static vector<CGene*> * getAllGenes() ;

	//GETTERS / SETTERS
	CVariant* getVariant(int position, string ref_allele, string novel_allele);
	CVariant* addVariant(CVariant*) throw (EGene);
	CVariant* remove(CVariant*);

	inline const vector<CVariant*> * getAllVariants() const { return _variants; }
	inline const string getName () const { return _name; };
	int getNumberOfProbandsWithVariants() const;

};

#endif /* CGENE_H_ */
