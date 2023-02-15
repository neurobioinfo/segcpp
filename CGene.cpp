/*
 * CGene.cpp
 *
 *  Created on: 2011-01-28
 *  Author: Alexandre Dionne-Laporte
 */

#include "CGene.h"
#include "CVariant.h"
#include "CSample.h"
#include "CChromosome.h"
#include "CVariantInSample.h"

#include <iostream>

using namespace std;

/** Global static pointer used to ensure a single instance of the class. */
map<string, CGene*> * CGene::_genes = NULL;

/** This function is called to create an instance of the class.
    Calling the constructor publicly is not allowed. The constructor
    is private and is only called by this Instance function.
 */
CGene* CGene::getSetGene(string name) {

	if (! _genes )   // Only allow one instance of class to be generated.
		_genes = new map<string, CGene*>();

	CGene* returnGene = NULL;
	map<string, CGene*>::iterator it = _genes->find(name);

	if ( it != _genes->end() ) {
		returnGene = it->second;
	}
	else {
		returnGene = new CGene (name);
		_genes->insert(pair<string, CGene*>(name,returnGene));
	}
	return returnGene;
}

vector<CGene*> * CGene::getAllGenes() {
	vector <CGene*> * genes = new vector <CGene*>();
	if (! _genes )   // Only allow one instance of class to be generated.
		_genes = new map<string, CGene*>();
	for ( map<string, CGene*>::iterator it=_genes->begin() ; it != _genes->end(); it++ ) {
		//cout << (*it).first << " => " << (*it).second << endl;
		genes->insert(genes->end(), (*it).second);
	}
	return genes;
}

CVariant* CGene::addVariant(CVariant* variant) throw (EGene) {
//	mtx_.lock();
	if ( getVariant(variant->getPosition(), variant->getRef_allele(), variant->getNovel_allele() ) == NULL ) {
		_variants->push_back( variant );
	}
//	else {
////		throw EGene ("Variant already present in gene", __FILE__, __LINE__);
//		return getVariant(variant->getPosition(), variant->getRef_allele(), variant->getNovel_allele() );
//	}
//	mtx_.unlock();
	return variant;
}

CVariant* CGene::remove(CVariant* variant) {
//	mtx_.lock();
	for (vector<CVariant*>::iterator it=_variants->begin(); it!=_variants->end(); ++it) {
		CVariant* compareToVariant = (*it);
		if (    (compareToVariant->getNovel_allele() == variant->getNovel_allele()) &&
				(compareToVariant->getRef_allele() == variant->getRef_allele()) &&
				(compareToVariant->getPosition() == variant->getPosition()) )
		{
			_variants->erase(it);
			break;
		}
	}
//	mtx_.unlock();
	return variant;
}

CVariant* CGene::getVariant(int position, string ref_allele, string novel_allele ) {
	CVariant* variant = NULL;
	for (vector<CVariant*>::const_iterator it=_variants->begin(); it!=_variants->end(); ++it) {
		CVariant* compareToVariant = (*it);
		if (    (compareToVariant->getNovel_allele() == novel_allele) &&
				(compareToVariant->getRef_allele() == ref_allele) &&
				(compareToVariant->getPosition() == position) )
		{
			variant = compareToVariant;
			break;
		}
	}
	return variant;
}

int CGene::getNumberOfProbandsWithVariants() const {

	if ( _variants == NULL ) { throw EGene ("Uninitialized variant list", __FILE__, __LINE__); }
	map <CSample*, int> sampleCounts;
	for ( vector<CVariant*>::const_iterator variant_it=_variants->begin() ; variant_it!=_variants->end() ; variant_it++ ) {
		CVariant *variant = *variant_it;
		const vector <CVariantInSample*> *variantsInSample = variant->getVariantInSamples();
		for ( vector <CVariantInSample*>::const_iterator variantInSample_it = variantsInSample->begin(); variantInSample_it != variantsInSample->end(); variantInSample_it++ ) {
			CVariantInSample* variantInSample = *variantInSample_it;
			CSample* sample = variantInSample->getSample();

			if ( ( (  sample->isProband() /* sample is proband and */
					  && ( variantInSample->getGenotype() == HOMOZYGOUS
						   || variantInSample->getGenotype() == HETEROZYGOUS
						   || variantInSample->getGenotype() == TRIZYGOUS
					  ) /* is a variant we want to show */
			       )
			       || CRunSettings::getInstance()->get_show_all_variants() /* or we show all variants regardless the genotype */
			     )
			   ) {
				sampleCounts[sample]++;
			}
		}
	}
	return sampleCounts.size();
}

std::ostream & operator << (std::ostream &stream , const CGene & gene) {
	return stream << "gene name: " << gene._name;
}
