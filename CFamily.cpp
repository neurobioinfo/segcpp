/*
 * CFamily.cpp
 *
 *  Created on: 2011-02-01
 *  Author: Alexandre Dionne-Laporte
 */

#include "CFamily.h"
#include "CSample.h"

using namespace std;

/** Global static pointer used to ensure a single instance of the class. */
map<string, CFamily*> * CFamily::_families = NULL;

/** This function is called to create an instance of the class.
    Calling the constructor publicly is not allowed. The constructor
    is private and is only called by this Instance function.
 */
CFamily* CFamily::getSetFamily(string id) {

	if (! _families ) {
		// Only allow one instance of class to be generated.
		_families = new map<string, CFamily*>();
	}

	CFamily* returnFamily = NULL;
	map<string, CFamily*>::iterator it = _families->find(id);

	if ( it != _families->end() ) {
		returnFamily = it->second;
	}
	else {
		returnFamily = new CFamily (id);
		_families->insert(pair<string,CFamily*>(id,returnFamily));
	}
	return returnFamily;
}

vector <CFamily*> * CFamily::getAllFamilies() {
	vector <CFamily*> * families = new vector <CFamily*>();
	for ( map<string, CFamily*>::iterator it=_families->begin() ; it != _families->end(); it++ ) {
		families->insert(families->end(), (*it).second);
	}
	return families;
}

int CFamily::countProbands () {
	int count = 0;
	for ( vector<CSample*>::iterator it=_samples->begin() ; it < _samples->end(); it++ ) {
		if ( (*it)->isProband() ) { count++; }
	}
	return count;
}
int CFamily::countMembers () {
	return _samples->size();
}

std::ostream & operator << (std::ostream &stream , const CFamily & family) {
	return stream << family._id;
}
