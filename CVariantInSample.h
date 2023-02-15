/*
 * CVariantInSample.h
 *
 *  Created on: 2011-02-15
 *  Author: Alexandre Dionne-Laporte
 */

#ifndef CVARIANTINSAMPLE_H_
#define CVARIANTINSAMPLE_H_

class CVariant;
class CSample;

#include "ITypes.h"

#include <string>
#include <map>

using namespace std;

class CVariantInSample {

protected:
	CSample * _sample;		    //snv_id
	CVariant * _variant; 		//variant
	double _genotype_quality;	//Genotype Quality (mpileup)
	VariantZygosity _genotype;  //Genotype zygosity
	int _depth;                 //variant call depth (DP field)
	int _genotyper_reads_ref;	//variant caller reads in ref (AD field)
	int _genotyper_reads_var;	//variant caller reads in var (AD field)
	int _total_reads_ref_plus;	//total reads in ref+ (ADS field)
	int _total_reads_ref_minus;	//total reads in ref- (ADS field)
	int _total_reads_var_plus;	//total reads in var+ (ADS field)
	int _total_reads_var_minus;	//total reads in var- (ADS field)
	map<string, string >* _formats; //Parsed required FORMAT fields

public:
	//CONSTRUCTORS  / DESTRUCTORS
	CVariantInSample(
			CSample * sample,
			CVariant * variant,
			double genotype_quality,
			VariantZygosity zygosity,
			int depth,
			int genotyper_reads_ref,
			int genotyper_reads_var,
			int total_reads_ref_plus,
			int total_reads_ref_minus,
			int total_reads_var_plus,
			int total_reads_var_minus
			);


	~CVariantInSample();

	//GETTERS
	inline CSample *getSample() const { return _sample;	}
	inline CVariant *getVariant() const { return _variant; }

	inline double getGenotypeQuality() const { return _genotype_quality; }
	inline VariantZygosity getGenotype () const { return _genotype; }

	/* aka DP filed */
	inline int getDepth () const { return _depth;}

	/* aka AD field */
	inline int getGenotyperReadsRef () const { return _genotyper_reads_ref;}
	inline int getGenotyperReadsVar () const { return _genotyper_reads_var;}
	inline int getGenotyperReads () const { return getGenotyperReadsRef() + getGenotyperReadsVar(); }
	inline double getGenotyperMutationFreq () const { return ((double)(getGenotyperReadsVar())) / ((double)(getGenotyperReads())); }

	/* aka ADS field */
	inline int getTotalReadsRefMinus() const {return _total_reads_ref_minus;}
	inline int getTotalReadsRefPlus() const {return _total_reads_ref_plus;}
	inline int getTotalReadsRef() const { return getTotalReadsRefPlus() + getTotalReadsRefMinus(); }
	inline int getTotalReadsVarMinus() const {return _total_reads_var_minus;}
	inline int getTotalReadsVarPlus() const {return _total_reads_var_plus;}
	inline int getTotalReadsVar() const {return getTotalReadsVarPlus() + getTotalReadsVarMinus(); }
	inline int getTotalReads() const {return getTotalReadsRef() + getTotalReadsVar(); }
	inline double getTotalMutationFreq() const { return ((double)(getTotalReadsVar())) / ((double)(getTotalReads())); }

	/* extra format fields */
	inline string setFormat (string key, string value) { (*_formats)[key]=value; return value;}
    inline const string getFormat (string key) const { return (*_formats)[key]; }
    inline const map<string, string >* getAllFormats () const { return _formats; }
};

#endif /* CVARIANTINSAMPLE_H_ */
