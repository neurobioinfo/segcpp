/*
 * CVCFFile.cpp
 *
 *  Created on: 2011-11-09
 *  Author: Alexandre Dionne-Laporte
 *  Inspired by vcflib, a simple C++ library for parsing and manipulating VCF files.
 */

#include "CVCFFile.h"
#include "CVariant.h"
#include "CSample.h"
#include "CFamily.h"
#include "UStringUtils.h"
#include "CChromosome.h"
#include "CGene.h"
#include "CVariantInSample.h"
#include "CRunSettings.h"
#include "ITypes.h"

#include <string>
#include <vector>
#include <iostream>
#include <float.h>

using namespace std;

const char* CVCFFile::ALLELE_FREQ = "AF";
const char* CVCFFile::DELS = "Dels";
const char* CVCFFile::STRAND_BIAS = "FS";
const char* CVCFFile::MAPPING_QUALITY = "MQ";
const char* CVCFFile::QUALITY_DEPTH = "QD";
const char* CVCFFile::CULPRIT = "culprit";

const char* CVCFFile::GENOTYPE = "GT";
const char* CVCFFile::ALLELE_DEPTH = "AD";
const char* CVCFFile::ALLELE_DEPTH_STRAND = "ADS";
const char* CVCFFile::DEPTH = "DP";
const char* CVCFFile::GENOTYPE_QUALITY = "GQ";

const char* CVCFFile::VARIANT_CLASS = "VC";
const char* CVCFFile::VARIANT_FUNCTION_TYPE = "VFT";
const char* CVCFFile::GENE = "Gene";
const char* CVCFFile::VARIANT_TYPE = "VT";
const char* CVCFFile::DETAILED_ANNOTATION = "DA";
//const char* CVCFFile::C_DBSNP_RS = "CdbSNP_rs";
//const char* CVCFFile::NC_DBSNP_RS = "NCdbSNP_rs";
//const char* CVCFFile::THOUSAND = "1KG";
//const char* CVCFFile::CG ="CG";
//const char* CVCFFile::EVS="ESP6500";
//const char* CVCFFile::SIFT = "SIFT";
//const char* CVCFFile::PPV2 = "PPv2";
//const char* CVCFFile::LRT = "LRT";
//const char* CVCFFile::MT = "MT";
//const char* CVCFFile::PHYLOP = "PhyloP";
//const char* CVCFFile::GERP = "GERP++";

void CVCFFile::openFile( string filename ) {
	//Assing variables
	_filename = filename;
	//Start parsing
	_file.open(_filename.c_str(), ifstream::in);
	_parseHeader();
}

vector<CVariant*>* CVCFFile::parseFile (
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
) throw (EVCFFile) {

	_min_proband_variant_reads = min_proband_variant_reads;
	_min_proband_coverage = min_proband_coverage;
	_min_proband_mutfreq = min_proband_mutfreq;
	_min_proband_genotype_quality = min_proband_genotype_quality;
	_min_familial_control_variant_reads = min_familial_control_variant_reads;
	_min_familial_control_coverage = min_familial_control_coverage;
	_min_familial_control_mutfreq = min_familial_control_mutfreq;
	_min_familial_control_genotype_quality = min_familial_control_genotype_quality;
	_min_external_control_variant_reads = min_external_control_variant_reads;
	_min_external_control_coverage = min_external_control_coverage;
	_min_external_control_mutfreq = min_external_control_mutfreq;
	_min_external_control_genotype_quality = min_external_control_genotype_quality;

	/* verbose option */
	int verbose_level = CRunSettings::getInstance()->get_verbose_level();
	if ( verbose_level > CRunSettings::VERBOSE_SILENT ) {
		cout << endl << output_seperator << endl;
		cout << "PARSING VCF FILE " << endl << "file name: " << _filename << endl <<
				"proband filters: " << _min_proband_variant_reads << "," << _min_proband_coverage << "," << _min_proband_mutfreq << "," << _min_proband_genotype_quality << " | " <<
				"fam controls filters: " << _min_familial_control_variant_reads << "," << _min_familial_control_coverage << "," << _min_familial_control_mutfreq << "," << _min_familial_control_genotype_quality << " | " <<
				"non fam controls filters: " << _min_external_control_variant_reads << "," << _min_external_control_coverage << "," << _min_external_control_mutfreq << "," << _min_external_control_genotype_quality << endl <<
				"[Pass Filter '.' | Total '-']" << endl << endl;
		cout.flush();
	}
	vector<CVariant*>* variants = new vector <CVariant*> ();
	vector<CVariant*>* newVariants = _getNextVariants();
	unsigned int filtered_total=0;
	while (!_finished) {
		unsigned int filtered_in_line=0;
		for ( vector<CVariant*>::iterator i = newVariants->begin() ; i != newVariants->end(); i++) {
			/* Adding previous line's variants*/
			if ( (*i)->getVariantInSamples()->empty() ) { delete *i;}

			else {
				variants->push_back(*i);

				/* checking for filtered variants */
				const vector<CVariantInSample*> *variantsInSample = (*i)->getVariantInSamples();
				for ( vector<CVariantInSample*>::const_iterator variantInSampleIt = variantsInSample->begin(); variantInSampleIt!=variantsInSample->end(); variantInSampleIt++ ) {
					if ( (*variantInSampleIt)->getGenotype() < 0 ) {
						filtered_in_line++;
					}
				}
			}
		}
		if ( filtered_in_line > newVariants->size()/2 ) {
			filtered_total++;
		}

		/* verbose option */
		if ( verbose_level > CRunSettings::VERBOSE_SILENT && _line_number % 100 == 0) {
			if ( newVariants->empty() ) {
				cout << "x" ;
			}
			else if ( filtered_in_line > newVariants->size()/2 ) {
				cout << "-" ;
			}
			else {
				cout << "." ;
			}

			if (_line_number % 10000 == 0) {
				cout << "    [" << variants->size()-filtered_total << " / " << variants->size() << "]" << endl;
			}
			cout.flush();
		}
		/* fetching next line's variants */
		newVariants = _getNextVariants();
	}
	/* verbose option */
	if ( verbose_level > CRunSettings::VERBOSE_SILENT ) {
		for ( int _line_number_copy = _line_number; _line_number_copy % 10000 != 0; _line_number_copy++) {
			if ( _line_number_copy % 100 == 0 ) {
				cout << " ";
			}
		}
		cout << "     [" << variants->size()-filtered_total << " / " << variants->size() << "]" << endl;
	}
	return variants;
}

void CVCFFile::_parseHeader() throw (EVCFFile) {
	string headerStr = "";

	while (std::getline(_file, _line)) {

		if (_line.substr(0,1) == "#") {
			headerStr += _line + '\n';

			string headerLine = _line;
			if (headerLine.substr(0,2) == "##") {
				// meta-information headerLines
				// TODO parse into map from info/format key to type
				// ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
				// ##FORMAT=<ID=CB,Number=1,Type=String,Description="Called by S(Sanger), M(UMich), B(BI)">
				size_t found = headerLine.find_first_of("=");
				string entryType = headerLine.substr(2, found - 2);
				// handle reference here, no "<" and ">" given
				//} else if (entryType == "reference") {
				size_t dataStart = headerLine.find_first_of("<");
				size_t dataEnd = headerLine.find_first_of(">");
				if (dataStart != string::npos && dataEnd != string::npos) {
					string entryData = headerLine.substr(dataStart + 1, dataEnd - dataStart - 1);

					// XXX bad; this will break if anyone ever moves the order
					// of the fields around to include a "long form" string
					// including either a = or , in the first or second field
					if (entryType == "INFO" || entryType == "FORMAT") {
						vector<string> fields = split(entryData, "=,");

						if (fields[0] != "ID") {
							throw EVCFFile(string("header parse error at:\n").
									append(fields[0]).append("!= \"ID\"\n").
									append(headerLine), __FILE__, __LINE__);
						}
						string id = fields[1];

						if (fields[2] != "Number") {
							throw EVCFFile(string("header parse error at:\n").
									append(fields[2]).append("!= \"Number\"\n").
									append(headerLine), __FILE__, __LINE__);
						}

						int variantFieldNumber;
						string numberstr = fields[3].c_str();
						variantFieldNumber = str2VariantFieldNumber (numberstr);

						if (fields[4] != "Type") {
							throw EVCFFile(string("header parse error at:\n").
									append(fields[4]).append("!= \"Type\"\n").
									append(headerLine), __FILE__, __LINE__);
						}
						VariantFieldType variantFieldType = str2VariantFieldType(fields[5]);

						if (fields[6] != "Description") {
							throw EVCFFile(string("header parse error at:\n").
									append(fields[6]).append("!= \"Description\"\n").
									append(headerLine), __FILE__, __LINE__);
						}
						string description;
						for ( unsigned int description_index=7 ; description_index<fields.size() ; description_index++ ) {
							description.append(fields[description_index]).append(",");
						}
						description = description.substr(0, description.size()-1);

						/* Inserting data */
						if (entryType == "INFO") {
							_infoCounts[id] = variantFieldNumber;
							_infoTypes[id] = variantFieldType;
							_infoDescriptions[id] = description;
						}
						else if (entryType == "FORMAT") {
							_formatCounts[id] = variantFieldNumber;
							_formatTypes[id] = variantFieldType;
							_formatDescriptions[id] = description;
						}
					}
				}
			}
			else if (headerLine.substr(0,1) == "#") {
				// field name headerLine
				vector<string> fields = split(headerLine, '\t');
				if (fields.size() > 8) {
					_sampleNames.resize(fields.size() - 9);
					copy(fields.begin() + 9, fields.end(), _sampleNames.begin());
				}
			}
		}
		else {
			// done with header
			if (headerStr.empty()) {
				throw EVCFFile(string("error: no VCFalt.size() header") , __FILE__, __LINE__);
			}
			_firstRecord = true;
			_finished = false;
			_line_number = 0;
			break;
		}
	}

	//Verifying if found all needed entries
	if ( _formatTypes.find(ALLELE_DEPTH_STRAND) == _formatTypes.end() && !CRunSettings::getInstance()->get_use_genotyper_reads()) {
		string errorSTR=string("error: no ").append(ALLELE_DEPTH_STRAND).append(" in FORMAT fields of file ");
		std::stringstream ss;
		ss << _filename.c_str();
		errorSTR.append(ss.str());
		throw EVCFFile(errorSTR , __FILE__, __LINE__);
	}
}

/**
 * Reads next line in VCF file and creates variants (one per alt allele)
 */
vector<CVariant*>* CVCFFile::_getNextVariants() throw (EVCFFile) {
	vector<CVariant*> * currentVariants = new vector<CVariant*>();
	if (_firstRecord) {
		/* then record was read from the parse method */
		_firstRecord = false;
	}
	else {
		std::getline(_file, _line);
	}

	if ( !_file.eof() ) {
		_line_number++;
		/*
		 * Line format:
		 *
		 * #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT [SAMPLE1 .. SAMPLEN]
		 */

		vector<string> fields = split(_line, '\t');
		string sequenceName = fields.at(0);
		unsigned long position;
		convert(fields.at(1).c_str(), position);
		string id = fields.at(2);
		string ref = fields.at(3);
		vector<string> alt = split(fields.at(4), ","); // a comma-separated list of alternate alleles
		double quality;
		convert(fields.at(5), quality);
		string filter = fields.at(6);

		/* Create variant(s) from this line; For them to exist to use in the samples */
		for ( unsigned int altNum=0; altNum < alt.size(); altNum++ ) {
			CVariant* variantTmp = CChromosome::getSetChromosome(sequenceName)->getVariant(position,ref,alt[altNum]);
			if ( variantTmp == NULL ) {
				VariantType variantType = _getVariantType(ref, alt[altNum]);
				variantTmp = new CVariant(
						CChromosome::getSetChromosome(sequenceName),
						position,
						ref,
						alt[altNum],
						variantType,
						filter
				);
			}
			currentVariants->push_back(variantTmp);
		}
		/* Parsing FORMAT fields */
		if (fields.size() > 8) {
			_parseFORMATFields(&fields, currentVariants );
		}
		else {
			throw new EVCFFile("Error with vcf file format",__FILE__, __LINE__);
		}

		/* Parsing INFO fields */
		_parseINFOFields(fields.at(7), currentVariants, alt.size());

		/* Add gene annotation */
		_addGeneAnnot(currentVariants);
	}
	else {
		_finished = true;
	}
	return currentVariants;
}

/**
 * Parse the variant INFO fields
 */
void CVCFFile::_parseINFOFields(string infoStr, vector<CVariant*> * currentVariants, unsigned int alt_num) throw (EVCFFile) {

	map<string, vector<string> > info_map;  // vector<string> allows for lists by Genotypes or Alternates
	vector<string> infofields = split(infoStr, ';');
	for (vector<string>::iterator info_it = infofields.begin(); info_it != infofields.end(); ++info_it) {
		if (*info_it == ".") {
			continue;
		}
		vector<string> info_hash = split(*info_it, '=');


		/*
		 * check if 1. field is required in output
		 */
		if (CRunSettings::getInstance()->get_is_required_info_fields(info_hash[0])) {

			if (info_hash.size() == 2) {
				info_map[info_hash[0]] = split(info_hash[1], ',');
				for ( unsigned int altNum=0; altNum < alt_num; altNum++ ) {
					CVariant* altVariant = currentVariants->at(altNum);
					
					/*
					 * check if 2. is present in at least one proband variant
					 */
					if ( altVariant->getNumberOfProbandsWithAltVariant() > 0 || CRunSettings::getInstance()->get_show_all_variants() ) {
						if ( getInfoCounts()[info_hash[0]] == ALLELE_NUMBER ) {
							vector <string> values = splitWithEmpty(info_hash[1], ",");
							if ( values.size() < altNum+1 ) {
								// TODO: throw exeption here
								cerr << endl << "Field in not concordant with file spec (allele #): (chr pos field): " <<
										altVariant->getChromosome()->getName() << " " <<
										altVariant->getPosition() << " " <<
										info_hash[0] << endl;
								altVariant->setInfo(info_hash[0], info_hash[1]);
							}
							else {
								if ( altVariant->getInfo(info_hash[0]).empty() ) {
									altVariant->setInfo(info_hash[0], values[altNum]);
								}
								else {
									//Adding value from other vcf to info field
									vector<string> field_values = splitStrDelim( altVariant->getInfo(info_hash[0]), "||");
									bool found=false;
									for ( vector<string>::iterator field_values_iter=field_values.begin(); field_values_iter!=field_values.end(); field_values_iter++ ) {
										if ( (*field_values_iter).compare(values[altNum]) == 0 ) {
											found=true;
											break;
										}
									}
									if ( !found ) {
										altVariant->setInfo(info_hash[0], altVariant->getInfo(info_hash[0]) + "||" + values[altNum]);
									}
								}
							}
						}
						else {
							if ( altVariant->getInfo(info_hash[0]).empty() ) {
								altVariant->setInfo(info_hash[0], info_hash[1]);
							}
							else {
								//Adding value from other vcf to info field
								vector<string> field_values = splitStrDelim( altVariant->getInfo(info_hash[0]), "||");
								bool found=false;
								for ( vector<string>::iterator field_values_iter=field_values.begin(); field_values_iter!=field_values.end(); field_values_iter++ ) {
									if ( (*field_values_iter).compare(info_hash[1]) == 0 ) {
										found=true;
										break;
									}
								}
								if ( !found ) {
									altVariant->setInfo(info_hash[0], altVariant->getInfo(info_hash[0]) + "||" + info_hash[1]);
								}
							}
						}
					}
					else {
//						cerr << *altVariant << "has no variant in proband" << endl;
					}
				}
			}
			else if (info_hash.size() == 1 && getInfoCounts()[info_hash[0]] == 0 ) {
				// Flag entry
				for ( unsigned int altNum=0; altNum < alt_num; altNum++ ) {
					CVariant* altVariant = currentVariants->at(altNum);
					altVariant->setInfo(info_hash[0], "Yes");
				}

			}
			else {
				vector<string> vec;
				vec.insert(vec.begin(), info_hash[0]);
				info_map[info_hash[0]] = vec;
			}
		}
	}
}

/** Parse the variant FORMAT fields */
void CVCFFile::_parseFORMATFields( vector<string> * fields,	vector<CVariant*> * currentVariants) throw (EVCFFile) {

	/* check if we have samples specified */

	vector<string> format = split(fields->at(8), ':');


	map<string, map<string, vector<string> > > samplesVariantsAttributes;  // vector<string> allows for lists by Genotypes or Alternates

	vector<string>::iterator sampleName = _sampleNames.begin();
	vector<string>::iterator sample = fields->begin() + 9;
	for (; sample != fields->end() && sampleName != _sampleNames.end(); ++sample, ++sampleName) {
		string& name = *sampleName;

//		if (*sample == ".") {
//			samplesVariantsAttributes.erase(name);
//			continue;
//		}

		vector<string> samplefields = split(*sample, ':');

/*		if ( *samplefields.begin() == "./.") {
			samplesVariantsAttributes[name][GENOTYPE].push_back("./.");
		}
		else */
		if ( ( *samplefields.begin() != "./.") && (samplefields.size() != format.size()) && (!CRunSettings::getInstance()->get_format_low_stringency()) ) {
			//We dont have the same number of fields in description and in sample's values
			string variantStr=string("samplefields ").append(*sample).append("   , format ").append(fields->at(8));
			throw EVCFFile(string("Error: Not the same number of fields in variant description and in sample's format values\n").append(variantStr), __FILE__, __LINE__);
		}
		else {
			/* Assinging FORMAT to the sample */
			vector<string>::iterator i = samplefields.begin();
			vector<string>::iterator f = format.begin();

			for (; f != format.end() && i != samplefields.end(); ++f, ++i) {
//				if (*i != ".") { samplesVariantsAttributes[name][*f] = split(*i, ','); }
//				to remove a bug and keep the "." genotype
				samplesVariantsAttributes[name][*f] = split(*i, ',');							
			}

		}
//		if (samplefields.size() != format.size()) {
//			/* allow "./." as a genotype, as many callers do it despite going against the sample specs. */
//			if (/*samplefields.size() == 1 &&*/ *samplefields.begin() == "./.")
//				samplesVariantsAttributes[name][GENOTYPE].push_back("./.");
//		}
//		else {
//			/* Assinging FORMAT to the sample */
//			vector<string>::iterator i = samplefields.begin();
//			vector<string>::iterator f = format.begin();
//			for (; f != format.end(); ++f, ++i) {
//				cout << "f " << *f << " i " << *i << endl;
//				samplesVariantsAttributes[name][*f] = split(*i, ',');
//			}
//		}
	}
	/* Assigning values to sample's corresponding Variant object */
	for (vector<string>::iterator sampleName = _sampleNames.begin(); sampleName != _sampleNames.end(); ++sampleName) {
		string& name = *sampleName;
//		cout << "name " << name << endl;
		if ( samplesVariantsAttributes.find(name) == samplesVariantsAttributes.end() ) {
			continue;
		}

		string genotype = samplesVariantsAttributes[name][GENOTYPE][0]; /* there is always only ONE genotype per sample */

//	 cout << "genotype: " << genotype << endl;

		CSample* sample = CSample::getSample(name); /* get sample object */
		if ( sample == NULL ) {
			/* This one is not defined in the ped file */
			continue;
		}
		else if ( (genotype.compare("./.") == 0) || (genotype.compare(".") == 0) ) {
			/**
			 * Todo add the variant as unknown with _addVariantInSample
			 */
			continue;
		}
		else {

			/* assing alleles from genotypes */
//			unsigned int allele1 = atoi(genotype.substr(0, genotype.find_first_of('/') ).c_str());
//			unsigned int allele2 = atoi(genotype.substr(genotype.find_first_of('/')+1, genotype.length()).c_str());
			int allele1 = -1;
			int allele2 = -1;
 			if ( genotype.find_first_of('/') ==  string::npos && genotype.find_first_of('|') !=  string::npos ) {
 				int allele1TMP = atoi(genotype.substr(0, genotype.find_first_of('|') ).c_str());
 				int allele2TMP = atoi(genotype.substr(genotype.find_first_of('|')+1, genotype.length()).c_str());
 				if  (allele1TMP > allele2TMP ) { allele1=allele2TMP; allele2=allele1TMP; }
 				else { allele1 = allele1TMP; allele2 = allele2TMP; }
 			}
 			else if (genotype.find_first_of('|') ==  string::npos && genotype.find_first_of('/') !=  string::npos ) {
                                int allele1TMP = atoi(genotype.substr(0, genotype.find_first_of('/') ).c_str());
                                int allele2TMP = atoi(genotype.substr(genotype.find_first_of('/')+1, genotype.length()).c_str());		
                                if  (allele1TMP > allele2TMP ) { allele1=allele2TMP; allele2=allele1TMP; }
                                else { allele1 = allele1TMP; allele2 = allele2TMP; }
			}
            /* DAN: special case: sex chromosomes in Illumina joint-genotyped data - male genotypes are just a single allele */
            else if ( genotype.find_first_of('/') == string::npos && genotype.find_first_of('|') == string::npos ) {
                        allele1 = atoi(genotype.c_str());
                        allele2 = allele1;
            }        
// 			cout << "alleles: " << allele1 << " : " << allele2 << endl;

			/* Loop over all variant alleles in the variant */
			if ( allele1 != allele2 && allele1 != 0) {
				currentVariants->at(allele1-1);
				currentVariants->at(allele2-1);
				samplesVariantsAttributes[name];
				/* first allele : trizygote */
				_addVariantInSample (sample, currentVariants->at(allele1-1), allele1, samplesVariantsAttributes[name], TRIZYGOUS);
				/* second allele */
				_addVariantInSample (sample, currentVariants->at(allele2-1), allele2, samplesVariantsAttributes[name], TRIZYGOUS);
			}
			else {

				/* One variant allele in this variant */
				VariantZygosity variantZygosity = ZYGOSITY_UNKNOWN;
				if (allele1==allele2) {
					if ( allele1 == 0 ) {
						variantZygosity = WILDTYPE;
					}
					else {
						variantZygosity = HOMOZYGOUS;
					}
				}
				else {
					variantZygosity = HETEROZYGOUS;
				}
				if ( variantZygosity == WILDTYPE ) {
					_addVariantInSample (sample, currentVariants->at(0), 0, samplesVariantsAttributes[name], variantZygosity);
				}
				else {
					_addVariantInSample (sample, currentVariants->at(allele2-1), allele2, samplesVariantsAttributes[name], variantZygosity);
				}
			}
		}

	}

	if (sampleName != _sampleNames.end()) {
		throw EVCFFile(string("error: more sample names in header than sample fields\n"), __FILE__, __LINE__);
	}
	if (sample != fields->end()) {
		throw EVCFFile(string("error: more sample fields than samples listed in header\n"), __FILE__, __LINE__);
	}
}

/** Add variant in sample */
CVariantInSample * CVCFFile::_addVariantInSample (
		CSample * sample,
		CVariant * variant,
		unsigned int var_allele_num,
		map<string, vector<string> > sampleVariantAttributes,
		VariantZygosity variantZygosity) throw (EVCFFile) {

//	cout << *sample << endl;
//	cout << variant->getPosition() << endl;



	CVariantInSample * variantInSample = NULL;

	/* Parsing DP field */
	int depth=0;
	if ( sampleVariantAttributes[DEPTH].size() > 0 ) {
		convert(sampleVariantAttributes[DEPTH][0], depth);
	}

	/* Parsing ADS field */
	int total_reads_ref_plus = 0;
	int total_reads_ref_minus = 0;
	int total_reads_ref = 0;
	int total_reads_var_plus = 0;
	int total_reads_var_minus = 0;
	int total_reads_var = 0;
	int total_reads = 0;
	double total_mutation_freq = 0.0;

//	cout << "stringency " << CRunSettings::getInstance()->get_format_low_stringency()  << endl;
//	cout << sampleVariantAttributes[ALLELE_DEPTH_STRAND].size() << endl;
	if ( sampleVariantAttributes.find(ALLELE_DEPTH_STRAND) != sampleVariantAttributes.end() && sampleVariantAttributes[ALLELE_DEPTH_STRAND][0].compare(".") != 0  ) {
		convert(sampleVariantAttributes[ALLELE_DEPTH_STRAND][0], total_reads_ref_plus);
		convert(sampleVariantAttributes[ALLELE_DEPTH_STRAND][1], total_reads_ref_minus);
		if ( variantZygosity == WILDTYPE ) {
			int total_reads_var_plus2 = 0;
			int total_reads_var_minus2 = 0;
			for ( unsigned int i=2; i < sampleVariantAttributes[ALLELE_DEPTH_STRAND].size(); i++) {
				convert(sampleVariantAttributes[ALLELE_DEPTH_STRAND][i++], total_reads_var_plus2);
				total_reads_var_plus += total_reads_var_plus2;
				convert(sampleVariantAttributes[ALLELE_DEPTH_STRAND][i], total_reads_var_minus2);
				total_reads_var_minus += total_reads_var_minus2;
			}
		}
		else {
			if ( sampleVariantAttributes[ALLELE_DEPTH_STRAND].size() == 4 ) {
				convert(sampleVariantAttributes[ALLELE_DEPTH_STRAND][2], total_reads_var_plus);
				convert(sampleVariantAttributes[ALLELE_DEPTH_STRAND][3], total_reads_var_minus);
			}
			else {
//				try {
					if ( sampleVariantAttributes[ALLELE_DEPTH_STRAND].size() > var_allele_num*2 ) {
					convert(sampleVariantAttributes[ALLELE_DEPTH_STRAND][var_allele_num*2], total_reads_var_plus);
					}
					else if ( !CRunSettings::getInstance()->get_format_low_stringency() ) {
						std::stringstream ss;
						ss << "error: cannot convert allele for variant " << *variant;
						throw EVCFFile(ss.str(), __FILE__, __LINE__);
					}

					if ( sampleVariantAttributes[ALLELE_DEPTH_STRAND].size() > (var_allele_num*2)+1) {
						convert(sampleVariantAttributes[ALLELE_DEPTH_STRAND][(var_allele_num*2)+1], total_reads_var_minus);
					}
					else if ( !CRunSettings::getInstance()->get_format_low_stringency() ) {
						std::stringstream ss;
						ss << "error: cannot convert allele for variant " << *variant;
						throw EVCFFile(ss.str(), __FILE__, __LINE__);
					}

			}
		}
		total_reads_ref = total_reads_ref_plus + total_reads_ref_minus;
		total_reads_var = total_reads_var_plus + total_reads_var_minus;
		total_reads = total_reads_ref + total_reads_var;
		total_mutation_freq = ((double)total_reads_var)/((double)total_reads);
	}
	else if ( !CRunSettings::getInstance()->get_use_genotyper_reads() && !CRunSettings::getInstance()->get_format_low_stringency() ) {
		string errMsg = string(ALLELE_DEPTH_STRAND).append(" field unavailable in sample ").append(sample->getName()).append(" pos ").append(variant->getChromosome()->getName()).append(" " + variant->getPosition() );
		throw EVCFFile(errMsg, __FILE__, __LINE__);
	}
	/* Parsing AD field */
	int genotyper_reads_ref = 0;
	int genotyper_reads_var = 0;
	int genotyper_reads = 0;
	double genotyper_mutation_freq = 0.0;
	if ( sampleVariantAttributes.find(ALLELE_DEPTH) != sampleVariantAttributes.end() && sampleVariantAttributes[ALLELE_DEPTH][0].compare(".") != 0 )
//             sampleVariantAttributes[ALLELE_DEPTH][0].compare(".") != 0 &&
//             sampleVariantAttributes[ALLELE_DEPTH].size() > var_allele_num ) 
        {

		convert(sampleVariantAttributes[ALLELE_DEPTH][0], genotyper_reads_ref);

		if ( variantZygosity == WILDTYPE ) {
			if ( sampleVariantAttributes[ALLELE_DEPTH].size() == 2 ) {
				convert(sampleVariantAttributes[ALLELE_DEPTH][1], genotyper_reads_var);
			}
		}
		else {
			if ( sampleVariantAttributes[ALLELE_DEPTH].size() == 2 ) {
				convert(sampleVariantAttributes[ALLELE_DEPTH][1], genotyper_reads_var);
			}
			else {
				if ( sampleVariantAttributes[ALLELE_DEPTH].size() > var_allele_num) {
					convert(sampleVariantAttributes[ALLELE_DEPTH][var_allele_num], genotyper_reads_var);
				}
				else if (!CRunSettings::getInstance()->get_format_low_stringency()) {
					std::stringstream ss;
					ss << ALLELE_DEPTH << " unavailable in sample " << sample->getName() << " @ " << variant->getChromosome()->getName() << ":" << variant->getPosition();
					throw EVCFFile(ss.str(), __FILE__, __LINE__);
				}
			}
		}
		genotyper_reads = genotyper_reads_ref + genotyper_reads_var;
		genotyper_mutation_freq = ((double)genotyper_reads_var)/((double)genotyper_reads);
	}
	else if (!CRunSettings::getInstance()->get_format_low_stringency()) {
		std::stringstream ss;
		ss << ALLELE_DEPTH << " unavailable in sample " << sample->getName() << " @ " << variant->getChromosome()->getName() << ":" << variant->getPosition();
		throw EVCFFile(ss.str(), __FILE__, __LINE__);
	}

	/**
	 * Filtering steps
	 */
	int compare_var_coverage = CRunSettings::getInstance()->get_use_genotyper_reads()? genotyper_reads_var:total_reads_var;
	int compare_coverage = CRunSettings::getInstance()->get_use_genotyper_reads()? genotyper_reads:total_reads;
	double compare_mutation_freq = CRunSettings::getInstance()->get_use_genotyper_reads()? genotyper_mutation_freq:total_mutation_freq;
	double genotype_qual = DBL_MIN;
	if ( sampleVariantAttributes[GENOTYPE_QUALITY].size() > 0 ) {
		convert(sampleVariantAttributes[GENOTYPE_QUALITY][0], genotype_qual);
	}
	/* to compare against */
	int sample_var_coverage_threshold = INT_MIN;
	int sample_coverage_threshold = INT_MIN;
	double sample_mutation_freq_threshold = DBL_MIN;
	double sample_genotype_qual_threshold = DBL_MIN;
	//First is proband sample
	if ( sample->isProband() ) {
		sample_coverage_threshold = _min_proband_coverage;
		sample_var_coverage_threshold = _min_proband_variant_reads;
		sample_mutation_freq_threshold = _min_proband_mutfreq;
		sample_genotype_qual_threshold = _min_proband_genotype_quality;
	}
	//controls
	else {
		//second is familial control
		if ( sample->getFamily()->countProbands() > 0) {
			sample_coverage_threshold = _min_familial_control_coverage;
			sample_var_coverage_threshold = _min_familial_control_variant_reads;
			sample_mutation_freq_threshold = _min_familial_control_mutfreq;
			sample_genotype_qual_threshold = _min_familial_control_genotype_quality;
		}
		//third is external control
		else {
			sample_coverage_threshold = _min_external_control_coverage;
			sample_var_coverage_threshold = _min_external_control_variant_reads;
			sample_mutation_freq_threshold = _min_external_control_mutfreq;
			sample_genotype_qual_threshold = _min_external_control_genotype_quality;
		}
	}

	/* wheter we need to check for filtered variant */
	if ( variant->isSNP() || CRunSettings::getInstance()->get_indel_filtering() ) {

		/* for Wildtype */
		if ( compare_coverage < sample_coverage_threshold      //check coverage
				|| (genotype_qual > 0 && genotype_qual < sample_genotype_qual_threshold) //check genotype quality
		) {
			variantZygosity = _filterZygocity(variantZygosity);
		}
		/* Variant allele */
		else if ( variantZygosity != WILDTYPE
				&& (  compare_mutation_freq < sample_mutation_freq_threshold  //check mutation frequency
						|| compare_var_coverage < sample_var_coverage_threshold //check number of variant reads
				)
		)
		{
			variantZygosity = _filterZygocity(variantZygosity);
		}
	}


	variantInSample = new CVariantInSample(
			sample,
			variant,
			genotype_qual,
			variantZygosity,
			depth,
			genotyper_reads_ref,
			genotyper_reads_var,
			total_reads_ref_plus,
			total_reads_ref_minus,
			total_reads_var_plus,
			total_reads_var_minus
	);


	const vector<string>* extra_format_fields = CRunSettings::getInstance()->get_extra_format_fields();
	for (vector<string>::const_iterator extra_format_fields_it = extra_format_fields->begin();
			extra_format_fields_it != extra_format_fields->end() ;
			extra_format_fields_it++ ) {
		//cout << "field : " << *extra_format_fields_it << endl;
		if ( sampleVariantAttributes[*extra_format_fields_it].empty() ) {
			continue;
		}
		else {
			// todo: check number elems in extra_format_fields_it
			//cout << "size of vector " << (sampleVariantAttributes[*extra_format_fields_it]).size();
			std::stringstream concat_format;

			for ( vector<string>::const_iterator format_it = sampleVariantAttributes[*extra_format_fields_it].begin();
					format_it != (sampleVariantAttributes[*extra_format_fields_it]).end();
					format_it++)
			{
				concat_format << (*format_it) << ",";
			}

			variantInSample->setFormat(*extra_format_fields_it, concat_format.str().substr (0,concat_format.str().size()-1));
		}
	}

//	for(map<string, string >::const_iterator it = variantInSample->getAllFormats()->begin(); it != variantInSample->getAllFormats()->end(); it++)
//	{
//	    string key = it->first;
//	    string value = it->second;
//	    //Do something
//	    //cout << key << ":" << value << "\t";
//	}
//	cout << endl;


//	if ( variantInSample->getGenotype() == WILDTYPE)
//	cout << variantInSample->getGenotyperReadsRef() << " " << variantInSample->getGenotyperReads() << endl;

	return variantInSample;
}

VariantZygosity CVCFFile::_filterZygocity ( VariantZygosity variantZygosity ) {
	if ( variantZygosity == HOMOZYGOUS ) {
		return HOMOZYGOUS_FILTERED;
	}
	else if ( variantZygosity == HETEROZYGOUS ) {
		return HETEROZYGOUS_FILTERED;
	}
	else if ( variantZygosity == WILDTYPE ) {
		return WILDTYPE_FILTERED;
	}
	else if ( variantZygosity == TRIZYGOUS ) {
		return TRIZYGOUS_FILTERED;
	}
	else {
		return ZYGOSITY_UNKNOWN;
	}
}

/** Add/Create gene associated to variant */
void CVCFFile::_addGeneAnnot(vector<CVariant*> * currentVariants) throw (EVCFFile) {
	for ( unsigned int altNum=0; altNum < currentVariants->size(); altNum++ ) {
		CVariant* altVariant = currentVariants->at(altNum);
		string geneName = altVariant->getInfo(CVCFFile::GENE);
		if ( !geneName.empty() ) {
			CGene *gene = CGene::getSetGene(geneName);
			altVariant->setGene(gene);
			gene->addVariant(altVariant);
		}
	}
}

/** Get Variant Type from vcf infos */
VariantType CVCFFile::_getVariantType (string ref, string var) {
	VariantType variantType = UNKNOWN;
	if ( var.compare("*") == 0) {
		variantType = INDEL;
	}        
	else if ( ref.size() == 1 && var.size() == 1 ) {
		variantType = SNV;
	}
	else if ( ref.size() != var.size() ) {
		variantType = INDEL;
	}
	return variantType;
}
