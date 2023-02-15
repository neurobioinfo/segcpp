/*
 * CVariant.h
 *
 *  Created on: 2011-01-28
 *  Author: Alexandre Dionne-Laporte
 */

#ifndef CVARIANT_H_
#define CVARIANT_H_

class CChromosome;
class CGene;
class CSample;
class CVariantInSample;

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <iostream>
#include <boost/thread/mutex.hpp>

#include "CGene.h"
#include "CVariantInSample.h"

using namespace std;

class CVariant {

private:
	boost::mutex mtx_;

protected:
	vector <CVariantInSample*> * _variantsInSamples;
	string _variant_caller;
	string _ref_allele; 	      //ref_allele
	string _novel_allele; 	      //muta
	CChromosome* _chromosome;     //chr
	int _position;			      //pos
	VariantType _type;		      //type
	string _filter;			      //Validity Tranches
	map<string, string >* _infos; //Parsed INFO fields
	CGene * _gene;                //The gene

public:

	//CONSTRUCTORS / DESCTRUCTORS
	CVariant (
		CChromosome* chromosome,
		int position,
		string ref_allele,
		string novel_allele,
		VariantType type,
		string filter
	);

	~CVariant();

	//GETTERS
    inline const vector <CVariantInSample*> * getVariantInSamples () const { return _variantsInSamples; }
    inline const CChromosome *getChromosome() const { return _chromosome; }
    inline const int getPosition() const { return _position; }
    inline const string getRef_allele() const { return _ref_allele; }
    inline const string getNovel_allele() const { return _novel_allele; }
    inline VariantType getType() const { return _type; }
    inline const string getFilter() const { return _filter; }
    inline const string getInfo (string key) const { return (*_infos)[key]; }
//    inline const map<string, string >* getInfos() const { return _infos; }
    inline const bool isSNP () const { return _ref_allele.length() == 1 && _novel_allele.length() == 1; }
    inline const CGene *getGene() const { return _gene; }

    //SETTERS
    void addVariantInSample (CVariantInSample* variantInSample);
    inline void setInfo ( string key, string value) {  (*_infos)[key] = value; }
    inline void setGene ( CGene* gene ) { _gene = gene;  }

    //Other functions implemented in cpp file
    const int getNumberOfProbandsWithAltVariant();
	friend std::ostream & operator << (std::ostream &, const CVariant &);

};

#endif /* CVARIANT_H_ */
