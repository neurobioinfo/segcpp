/*
 * CVCFFile.h
 *
 *  Created on: 2011-11-09
 *  Author: Alexandre Dionne-Laporte
 */

#ifndef CVCF_FILE_H_
#define CVCF_FILE_H_

#include "ITypes.h"
#include "CRunSettings.h"

class CVariant;
class CSample;
class CVariantInSample;

#include <fstream>
#include <map>
#include <vector>
#include <sstream>
#include <iostream>

using namespace std;

class EVCFFile : public std::exception {
private:
	string _msg;
	string _file;
	string _line;

public:
	EVCFFile(const string msg, const string file, const int line) : _msg(msg), _file(file) {
		std::stringstream ss;
		ss << line;
		_line = ss.str();
	};
	virtual ~EVCFFile() throw() {}
//	virtual const char* what() const throw() {
//		string ret;
//		ret.append(_msg).append("\nIn file ").append(_file).append(", line ").append(_line);
//		cerr << ret << endl;
//		return ret.c_str();
//	}
	const string message() {
		std::stringstream ss;
		ss << _msg << endl << "In file " << _file << ", line " << _line << endl;;
		return ss.str();
	}

};

class CVCFFile {
public:
	CVCFFile() {};
	~CVCFFile() {
		/* VERBOSE output */
//		if (CRunSettings::getInstance()->get_verbose_level() > CRunSettings::VERBOSE_LOG)
//			cerr << "CVCFFile::~CVCFFile" << endl;
	};

	/* Regular methods accessible outside the class */
	void openFile( string filename );
	vector<CVariant*>* parseFile (
			int min_proband_variant_reads,
			int min_proband_coverage,
			double min_proband_mutfreq,
			int min_proband_genotype_quality,
			int min_familial_control_variant_reads,
			int min_familial_control_coverage,
			double min_familial_control_mutfreq,
			int min_familial_control_genotype_quality,
			int min_external_control_variant_reads,
			int min_external_control_coverage,
			double min_external_control_mutfreq,
			int min_external_control_genotype_quality
	) throw (EVCFFile);

	/* Getters */
	inline map<string, VariantFieldType> getInfoTypes() { return _infoTypes; }
	inline map<string, int> getInfoCounts() {return _infoCounts; }
	inline map<string, string> getInfoDescriptions() { return _infoDescriptions; }
	inline map<string, VariantFieldType> getFormatTypes() { return _formatTypes; }
	inline map<string, int> getFormatCounts() { return _formatCounts; }
	inline map<string, string> getFormatDescriptions() { return _formatDescriptions; }
	inline vector<string> getSampleNames() { return _sampleNames; }
	inline int getLineNumber () {return _line_number; }

	/* CONSTANTS */
	static const char* ALLELE_FREQ;
	static const char* DELS;
	static const char* STRAND_BIAS;
	static const char* MAPPING_QUALITY;
	static const char* QUALITY_DEPTH;
	static const char* CULPRIT;

	static const char* GENOTYPE;
	static const char* ALLELE_DEPTH;
	static const char* ALLELE_DEPTH_STRAND;
	static const char* DEPTH;
	static const char* GENOTYPE_QUALITY;

	static const char* VARIANT_CLASS;
	static const char* VARIANT_FUNCTION_TYPE;
	static const char* GENE;
	static const char* VARIANT_TYPE;
	static const char* DETAILED_ANNOTATION;
	static const char* C_DBSNP_RS;
	static const char* NC_DBSNP_RS;
	static const char* THOUSAND;
	static const char* CG;
	static const char* EVS;
	static const char* SIFT;
	static const char* PPV2;
	static const char* LRT;
	static const char* MT;
	static const char* PHYLOP;
	static const char* GERP;
	static const char* MIRNA;
	static const char* MIRNATARGET;

private:
	/* Variables */
	ifstream _file;
	string _filename;
	string _line; // the current line
	int _line_number; // the current line number
	bool _parsedHeader;
	bool _firstRecord;
	bool _finished;

	map<string, VariantFieldType> _infoTypes;
	map<string, int> _infoCounts;
	map<string, string> _infoDescriptions;
	map<string, VariantFieldType> _formatTypes;
	map<string, int> _formatCounts;
	map<string, string> _formatDescriptions;
	vector<string> _sampleNames;

	/* Filtration variables specific to the vcf file */
	int _min_proband_variant_reads;
	int _min_familial_control_variant_reads;
	int _min_external_control_variant_reads;

	int _min_proband_coverage;
	int _min_familial_control_coverage;
	int _min_external_control_coverage;

	double _min_proband_mutfreq;
	double _min_familial_control_mutfreq;
	double _min_external_control_mutfreq;

	int _min_proband_genotype_quality;
	int _min_familial_control_genotype_quality;
	int _min_external_control_genotype_quality;

	/* Functions */
	void _parseHeader() throw (EVCFFile) ;
	vector<CVariant*>* _getNextVariants() throw (EVCFFile) ;
	void _parseINFOFields(string infoStr, vector<CVariant*> * currentVariants, unsigned int alt_num) throw (EVCFFile);
	void _parseFORMATFields(vector<string> * fields, vector<CVariant*> * currentVariants) throw (EVCFFile);
	void _addGeneAnnot(vector<CVariant*> * currentVariants) throw (EVCFFile);
	CVariantInSample * _addVariantInSample (CSample * sample, CVariant * variant, unsigned int var_allele_num, map<string, vector<string> > sampleVariantsAttributes, VariantZygosity variantZygosity) throw (EVCFFile);
	VariantZygosity _filterZygocity ( VariantZygosity variantZygosity );
	VariantType _getVariantType (string ref, string var);

};

#endif /* CVCFFILE_H_ */
