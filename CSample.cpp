/*
 * CSample.cpp
 *
 * Created on: 2011-02-01
 * Author: Alexandre Dionne-Laporte
 */

#include "CSample.h"
#include "CFamily.h"
#include "CRunSettings.h"

#include <iostream>

using namespace std;

/** Global static pointer used to ensure a single instance of the class. */
map<string, CSample*> * CSample::_samples = new map<string, CSample*>();

/** This function is called to create an instance of the class.
    Calling the constructor publicly is not allowed. The constructor
    is private and is only called by this Instance function.
 */
CSample* CSample::getSample(string name) {
	CSample* returnCSample = NULL;
	map<string, CSample*>::iterator it = _samples->find(name);

	if ( it != _samples->end() ) {
		returnCSample = it->second;
	}
	return returnCSample;
}

vector<CSample*> CSample::getSamples() {
	vector<CSample*> values;
	for(map<string, CSample*>::iterator it = _samples->begin(); it != _samples->end(); ++it) {
		values.push_back(it->second);
	}
	return values;
}


CSample::CSample(string name, bool isProband, CFamily * family) :
		_name(name), _isProband(isProband), _family(family) {
	_variantsInSample = new vector<CVariantInSample*>();
	(*_samples)[name]=this;
};

CSample::~CSample() {
	_variantsInSample->clear();
	delete _variantsInSample;
	_samples->erase((*_samples).find(_name)) ;
//	if (CRunSettings::getInstance()->get_verbose_level() > CRunSettings::VERBOSE_LOG)
//		cerr << "CSample::~CSample" << endl;
}

std::ostream & operator << (std::ostream &stream, const CSample & sample) {
	return stream << "name: " << sample._name << ", is proband: " << sample._isProband << ", family: " <<  *(sample._family);
}
