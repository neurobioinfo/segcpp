/*
 * CChromosome.cpp
 *
 *  Created on: 2011-01-28
 *  Author: Alexandre Dionne-Laporte
 */

#include "CChromosome.h"
#include "CVariant.h"

#include <iostream>
#include <boost/thread/mutex.hpp>

using namespace std;

/** Global static pointer used to ensure a single instance of the class. */
map<string, CChromosome*> * CChromosome::_chromosomes = NULL;

/** This function is called to create an instance of the class.
    Calling the constructor publicly is not allowed. The constructor
    is private and is only called by this Instance function.
 */
CChromosome* CChromosome::getSetChromosome(string name) {

	if (! _chromosomes )   // Only allow one instance of class to be generated.
		_chromosomes = new map<string, CChromosome*>();

	CChromosome* returnChromosome = NULL;
	map<string, CChromosome*>::iterator it = _chromosomes->find(name);

	if ( it != _chromosomes->end() ) {
		returnChromosome = it->second;
	}
	else {
		returnChromosome = new CChromosome (name);
		_chromosomes->insert(pair<string, CChromosome*>(name,returnChromosome));
	}
	return returnChromosome;
}

vector<CChromosome*> * CChromosome::getAllChromosomes() {
	vector <CChromosome*> * chromosomes = new vector <CChromosome*>();
	if (! _chromosomes )   // Only allow one instance of class to be generated.
		_chromosomes = new map<string, CChromosome*>();
	for ( map<string, CChromosome*>::iterator it=_chromosomes->begin() ; it != _chromosomes->end(); it++ ) {
		//cout << (*it).first << " => " << (*it).second << endl;
		chromosomes->insert(chromosomes->end(), (*it).second);
	}
	return chromosomes;
}

CVariant* CChromosome::addVariant(CVariant* variant) {
//	mtx_.lock();
	if ( getVariant(variant->getPosition(), variant->getRef_allele(), variant->getNovel_allele() ) == NULL ) {
		_variants->insert( pair<int,CVariant*>(variant->getPosition(), variant) );
	}
	else {
		throw EChromosome ("Variant already present in chromosome", __FILE__, __LINE__);
	}
//	mtx_.unlock();
	return variant;
}

CVariant* CChromosome::remove(CVariant* variant) {
//	mtx_.lock();
	pair<multimap<int,CVariant*>::iterator,multimap<int, CVariant*>::iterator> ret;
	multimap<int,CVariant*>::iterator it;
	ret = _variants->equal_range(variant->getPosition());
	for (it=ret.first; it!=ret.second; ++it) {
		CVariant* compareToVariant = (*it).second;
		if (    (compareToVariant->getNovel_allele() == variant->getNovel_allele()) &&
				(compareToVariant->getRef_allele() == variant->getRef_allele()) )
		{
			_variants->erase(it);
			break;
		}
	}
//	mtx_.unlock();
	return variant;
}

CVariant* CChromosome::getVariant(int position, string ref_allele, string novel_allele ) {

	CVariant* variant = NULL;
	pair<multimap<int,CVariant*>::iterator,multimap<int, CVariant*>::iterator> ret;
	multimap<int,CVariant*>::iterator it;
	ret = _variants->equal_range(position);
	for (it=ret.first; it!=ret.second; ++it) {
		CVariant* compareToVariant = (*it).second;
		if (    (compareToVariant->getNovel_allele() == novel_allele) &&
				(compareToVariant->getRef_allele() == ref_allele) )
		{
			variant = compareToVariant;
			break;
		}
	}
	return variant;
}

const vector<CVariant*> * CChromosome::getAllVariants() {
	vector <CVariant*> * variants = new vector <CVariant*>();
	for ( map<int, CVariant*>::iterator it=_variants->begin() ; it != _variants->end(); it++ ) {
		variants->insert(variants->end(), (*it).second);
	}
	return variants;
}

std::ostream & operator << (std::ostream &stream , const CChromosome & chromosome) {
	return stream << chromosome._name;
}
