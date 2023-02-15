/*
 * CFamily.h
 *
 *  Created on: 2011-02-01
 *  Author: Alexandre Dionne-Laporte
 */

#ifndef CFAMILY_H_
#define CFAMILY_H_

class CSample;
#include "CRunSettings.h"

#include <map>
#include <vector>
#include <string>
#include <iosfwd>
#include <iostream>

using namespace std;

class CFamily {

	friend std::ostream & operator << (std::ostream &, const CFamily &);

protected:
	string _id;
	vector<CSample*> * _samples;

private:
	static map<string, CFamily*> * _families;
	CFamily(string id) : _id(id) { _samples = new vector<CSample*>(); };
	~CFamily() {
		_samples->clear();
		delete _samples;
//		if (CRunSettings::getInstance()->get_verbose_level() > CRunSettings::VERBOSE_LOG)
//			cerr << "CFamily::~CFamily" << endl;
	}

public:
	static CFamily* getSetFamily(string);
	static vector<CFamily*> * getAllFamilies() ;
	inline void addSample (CSample* sample) { _samples->insert(_samples->end(), sample); };
	inline vector<CSample*> * getSamples () const { return _samples; };
	int countProbands () ;
	int countMembers ();

	//GETTERS
	inline string getID() { return _id; };

};

#endif /* CFAMILY_H_ */
