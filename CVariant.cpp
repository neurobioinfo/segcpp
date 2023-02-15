/*
 * CVariant.cpp
 *
 *  Created on: 2011-01-28
 *  Author: Alexandre Dionne-Laporte
 */

#include "CVariant.h"
#include "CVariantInSample.h"
#include "CChromosome.h"
#include "CGene.h"
#include "CSample.h"
#include "CFamily.h"
#include "CVCFFile.h"
#include "CRunSettings.h"
#include "UStringUtils.h"

#include <set>
#include <iostream>

using namespace std;

/**
 * Constructors / Desctructors
 */
CVariant::CVariant (
	CChromosome* chromosome,
	int position,
	string ref_allele,
	string novel_allele,
	VariantType type,
	string filter
	) :
	_ref_allele(ref_allele),
	_novel_allele(novel_allele),
	_chromosome(chromosome),
	_position(position),
	_type (type),
	_filter(filter)

{
	_infos = new map<string, string >();
	_variantsInSamples = new vector <CVariantInSample*>();
	_chromosome->addVariant(this);
	_gene = NULL;
}

CVariant::~CVariant() {
	_infos->clear();
	delete _infos;
	_variantsInSamples->clear();
	delete _variantsInSamples;
	_chromosome->remove(this);
	if ( _gene != NULL ) { _gene->remove(this); }
	/* VERBOSE output */
//	if (CRunSettings::getInstance()->get_verbose_level() > CRunSettings::VERBOSE_LOG)
//		cerr << "CVariant::~CVariant" << endl;
};

void CVariant::addVariantInSample (CVariantInSample* variantInSample) {
	_variantsInSamples->push_back( variantInSample );
}

const int CVariant::getNumberOfProbandsWithAltVariant() {
	int probandVariants=0;
	const vector<CVariantInSample*> *variantsInSample = getVariantInSamples();
	for ( vector<CVariantInSample*>::const_iterator variantInSampleIt = variantsInSample->begin(); variantInSampleIt!=variantsInSample->end(); variantInSampleIt++ ) {
		if (    (*variantInSampleIt)->getSample()->isProband() &&
				(
						(*variantInSampleIt)->getGenotype() == HETEROZYGOUS ||
						(*variantInSampleIt)->getGenotype() == HOMOZYGOUS ||
						(*variantInSampleIt)->getGenotype() == TRIZYGOUS ||
						(*variantInSampleIt)->getGenotype() == HETEROZYGOUS_FILTERED ||
						(*variantInSampleIt)->getGenotype() == HOMOZYGOUS_FILTERED ||
						(*variantInSampleIt)->getGenotype() == TRIZYGOUS_FILTERED
				)
		) {
			probandVariants++;
		}
	}
	return probandVariants;
}

std::ostream & operator << (std::ostream &stream , const CVariant & variant) {
	set<CFamily*> printedFamilies;
	const vector<CVariantInSample*> * variantInSample = variant.getVariantInSamples();

	//loop over each samples in the current variant
	for ( vector<CVariantInSample*>::const_iterator var_sample_it=variantInSample->begin() ; var_sample_it < variantInSample->end(); var_sample_it++ ) {
		CVariantInSample * variantInSample = *var_sample_it;
		CSample * sample = variantInSample->getSample();
		if ( ( (  sample->isProband() /* sample is proband and */
			      && ( variantInSample->getGenotype() == HOMOZYGOUS
			    	   || variantInSample->getGenotype() == HETEROZYGOUS
			    	   || variantInSample->getGenotype() == TRIZYGOUS
			         ) /* is a variant we want to show */
		       )
			 || CRunSettings::getInstance()->get_show_all_variants() /* or we show all variants regardless the genotype */
			 )
		   && (printedFamilies.find(sample->getFamily()) == printedFamilies.end()) /* and we havent print this variant from this family */
		   )
		{

			CSample * _sample = variantInSample->getSample();
			CVariant * _variant = variantInSample->getVariant();

			/* Prints variant shared information */
			stream << _sample->getFamily()->getID() << "\t";
			stream << _variant->getChromosome()->getName() << "\t";
			stream << _variant->getPosition() << "\t";
			stream << _variant->getRef_allele() << "\t";
			stream << _variant->getNovel_allele() << "\t";

//			stream << _variant->getInfo(CVCFFile::VARIANT_CLASS) << "\t";
//			stream << _variant->getInfo(CVCFFile::GENE) << "\t";
//			stream << _variant->getInfo(CVCFFile::VARIANT_FUNCTION_TYPE) << "\t";
//			stream << _variant->getInfo(CVCFFile::VARIANT_TYPE) << "\t";
//			stream << _variant->getInfo(CVCFFile::DETAILED_ANNOTATION) << "\t";

			if ( _variant->getGene() != NULL ) {
				stream << _variant->getGene()->getNumberOfProbandsWithVariants() << "\t";
			}
			else {
				stream << "\t";
			}

//			stream << _variant->getInfo(CVCFFile::C_DBSNP_RS) << "\t";
//			stream << _variant->getInfo(CVCFFile::NC_DBSNP_RS) << "\t";
//			stream << _variant->getInfo(CVCFFile::THOUSAND) << "\t";
//			stream << _variant->getInfo(CVCFFile::CG) << "\t";
//			stream << _variant->getInfo(CVCFFile::EVS) << "\t";
//			stream << _variant->getInfo(CVCFFile::SIFT) << "\t";
//			stream << _variant->getInfo(CVCFFile::PPV2) << "\t";
//			stream << _variant->getInfo(CVCFFile::LRT) << "\t";
//			stream << _variant->getInfo(CVCFFile::MT) << "\t";
//			stream << _variant->getInfo(CVCFFile::PHYLOP) << "\t";
//			stream << _variant->getInfo(CVCFFile::GERP) << "\t";

			const vector<string>* extra_info_fields = CRunSettings::getInstance()->get_extra_info_fields();
			for (vector<string>::const_iterator extra_info_fields_it = extra_info_fields->begin();
					extra_info_fields_it != extra_info_fields->end() ;
					extra_info_fields_it++ ) {
				if ( _variant->getInfo(*extra_info_fields_it).empty() ) {
					stream << "\t";
				}
				else {
					stream << _variant->getInfo(*extra_info_fields_it) << "\t";
				}
			}

//			stream << _variant->getInfo(CVCFFile::ALLELE_FREQ) << "\t";
//			stream << _variant->getInfo(CVCFFile::DELS) << "\t";
//			stream << _variant->getInfo(CVCFFile::STRAND_BIAS) << "\t";
//			stream << _variant->getInfo(CVCFFile::MAPPING_QUALITY) << "\t";
//			stream << _variant->getInfo(CVCFFile::QUALITY_DEPTH) << "\t";
			stream << _variant->getFilter() << "\t";
//			stream << _variant->getInfo(CVCFFile::CULPRIT) << "\t";

//			string external_control_names_w_variant;
			vector <string> external_control_names_w_variant;
			vector <string> proband_names_w_variant;
			vector <string> familial_control_names_w_variant;
			int total_subjects_filtered=0;
			int total_subjects_wildtype=0;
			int total_subjects_w_variant=0;

			int total_affected_w_variant=0;
			int total_controls_w_variant=0;

			int external_controls_filtered=0;
			int external_controls_wildtype=0;
			int external_controls_w_variant=0;
			int external_controls_w_variant_homo=0;

			int family_proband_filtered=0;
			int family_proband_wildtype=0;
			int family_proband_w_variant=0;

			int family_controls_filtered=0;
			int family_controls_wildtype=0;
			int family_controls_w_variant=0;

			int sum_total_alleles = 0;
			int sum_variant_alleles = 0;

			const vector <CVariantInSample*> * variantInSamples = _variant->getVariantInSamples();
			for ( vector <CVariantInSample*>::const_iterator var_sample_it=variantInSamples->begin() ; var_sample_it < variantInSamples->end(); var_sample_it++ ) {
				CVariantInSample * variantInSample = *var_sample_it;
				CSample * sample = variantInSample->getSample();

				//stats for Filtered calls
				if ( (variantInSample->getGenotype() == HOMOZYGOUS_FILTERED )
					 || (variantInSample->getGenotype() == HETEROZYGOUS_FILTERED )
					 || (variantInSample->getGenotype() == TRIZYGOUS_FILTERED )
					 || (variantInSample->getGenotype() == WILDTYPE_FILTERED )
				   )
				{
					total_subjects_filtered++;
					if ( sample->getFamily()->getID() == _sample->getFamily()->getID() ) {
						//for family members
						if ( sample->isProband() ) {
							family_proband_filtered++;
						}
						else {
							family_controls_filtered++;
						}
					}
					else {
						bool found_proband = false;
						vector<CSample*> * family_samples = sample->getFamily()->getSamples();
						for ( vector<CSample*>::iterator family_samples_it=family_samples->begin(); family_samples_it!=family_samples->end(); family_samples_it++) {
							if ( (*family_samples_it)->isProband() ) {
								found_proband = true;
								break;
							}
						}
						if (!found_proband) { external_controls_filtered++; }
					}
				}

				//stats for Wildtype calls, REMOVING FILTERED STUFF
				if ( variantInSample->getGenotype() == WILDTYPE )
				{
					total_subjects_wildtype++;
					if ( sample->getFamily()->getID() == _sample->getFamily()->getID() ) {
						//for family members
						if ( sample->isProband() ) {
							family_proband_wildtype++;
						}
						else {
							family_controls_wildtype++;
						}
					}
					else {
						bool found_proband = false;
						vector<CSample*> * family_samples = sample->getFamily()->getSamples();
						for ( vector<CSample*>::iterator family_samples_it=family_samples->begin(); family_samples_it!=family_samples->end(); family_samples_it++) {
							if ( (*family_samples_it)->isProband() ) {
								found_proband = true;
								break;
							}
						}
						if (!found_proband) { external_controls_wildtype++; }
					}
				}

				//Stats for variant calls
				if ( (variantInSample->getGenotype() == HOMOZYGOUS )
					 || (variantInSample->getGenotype() == HETEROZYGOUS )
					 || (variantInSample->getGenotype() == TRIZYGOUS )
				   )
				{
					total_subjects_w_variant++;

					if ( sample->isProband() ) {
						total_affected_w_variant++;
						if ( sample->getFamily()->getID() == _sample->getFamily()->getID() ) {
							family_proband_w_variant++;
							proband_names_w_variant.push_back(sample->getName());
						}

					}
					else {
						total_controls_w_variant++;

						if ( sample->getFamily()->getID() == _sample->getFamily()->getID() ) {
							family_controls_w_variant++;
							familial_control_names_w_variant.push_back(sample->getName());
						}
						else {
							//external controls
							bool found_proband = false;
							vector<CSample*> * family_samples = sample->getFamily()->getSamples();
							for ( vector<CSample*>::iterator family_samples_it=family_samples->begin(); family_samples_it!=family_samples->end(); family_samples_it++) {
								if ( (*family_samples_it)->isProband() ) {
									found_proband = true;
									break;
								}
							}
							if (!found_proband) {
								//we are in the case of an external control
								external_control_names_w_variant.push_back(sample->getName());
								external_controls_w_variant++;
								if ( variantInSample->getGenotype() == HOMOZYGOUS ) {
									external_controls_w_variant_homo++;
								}
							}
						}
					}
				}

				/* Counting Alleles */
				if (variantInSample->getGenotype() == WILDTYPE) {
					sum_total_alleles += 2;
					sum_variant_alleles += 0;
				}
				else if (variantInSample->getGenotype() == HETEROZYGOUS) {
					sum_total_alleles += 2;
					sum_variant_alleles += 1;
				}
				else if (variantInSample->getGenotype() == HOMOZYGOUS) {
					sum_total_alleles += 2;
					sum_variant_alleles += 2;
				}
				else if (variantInSample->getGenotype() == TRIZYGOUS) {
					sum_total_alleles += 2;
					sum_variant_alleles += 2;
				}

			}

			/* Printing allele frequencies */
			stream << sum_total_alleles << "\t";
			stream << (double)sum_variant_alleles / sum_total_alleles << "\t";

			/* Proband names */
			string proband_names_STR="";
			if ( proband_names_w_variant.size() > CRunSettings::getInstance()->get_show_probands_with_variants() ) {
				proband_names_STR = "...";
			}
			else {
				for ( unsigned int proband_name_index=0; proband_name_index<proband_names_w_variant.size();	proband_name_index++) {
					proband_names_STR.append(proband_names_w_variant.at(proband_name_index).append(","));
				}
				if ( proband_names_STR.length() > 0 ) {	proband_names_STR = proband_names_STR.substr(0, proband_names_STR.size()-1); }
			}
			stream << proband_names_STR << "\t";
			/* Familial control names */
			string familial_control_names_STR="";
			if ( familial_control_names_w_variant.size() > CRunSettings::getInstance()->get_show_familial_controls_with_variants() ) {
				familial_control_names_STR = "...";
			}
			else {
				for ( unsigned int familial_control_name_index=0; familial_control_name_index<familial_control_names_w_variant.size(); familial_control_name_index++) {
					familial_control_names_STR.append(familial_control_names_w_variant.at(familial_control_name_index).append(","));
				}
				if ( familial_control_names_STR.length() > 0 ) { familial_control_names_STR = familial_control_names_STR.substr(0, familial_control_names_STR.size()-1); }
			}
			stream << familial_control_names_STR << "\t";

			/* External control names */
			string external_control_names_STR="";
			if ( external_control_names_w_variant.size() > CRunSettings::getInstance()->get_show_external_controls_with_variants() ) {
				external_control_names_STR = "...";
			}
			else {
				for ( unsigned int control_name_index=0; control_name_index<external_control_names_w_variant.size(); control_name_index++) {
					external_control_names_STR.append(external_control_names_w_variant.at(control_name_index).append(","));
				}
				if ( external_control_names_STR.size() > 0 ) { external_control_names_STR = external_control_names_STR.substr(0, external_control_names_STR.size()-1); }
			}
			stream << external_control_names_STR << "\t";

			/* variant counts columns */
			int family_total_w_variant = family_proband_w_variant + family_controls_w_variant;
			stream << total_subjects_filtered << "\t" <<
					total_subjects_wildtype << "\t" <<
					total_subjects_w_variant << "\t" <<
					total_affected_w_variant << "\t" <<
					total_controls_w_variant << "\t" <<
					external_controls_filtered << "\t" <<
					external_controls_wildtype << "\t" <<
					external_controls_w_variant << "\t" <<
					external_controls_w_variant_homo << "\t" <<
					family_total_w_variant << "\t" <<
					family_proband_filtered << "\t" <<
					family_proband_wildtype << "\t" <<
					family_proband_w_variant << "\t" <<
					family_controls_filtered << "\t" <<
					family_controls_wildtype << "\t" <<
					family_controls_w_variant << "\t";

			string family_members_str;
			string family_members_info;

			for ( vector<CSample*>::iterator family_members_it=_sample->getFamily()->getSamples()->begin() ; family_members_it < _sample->getFamily()->getSamples()->end(); family_members_it++ ) {
				CSample* family_member = *family_members_it;
				bool found_individual = false;
				family_members_str.append(family_member->getName());
				family_members_str.append(",");

				for ( vector <CVariantInSample*>::const_iterator variant_in_samples_it = _variant->getVariantInSamples()->begin() ;
						variant_in_samples_it < _variant->getVariantInSamples()->end() ;
						variant_in_samples_it++ ) {
					CVariantInSample* variant_in_sample = *variant_in_samples_it;
					if ( variant_in_sample->getSample() == family_member ) {
						found_individual = true;

						double GQ = variant_in_sample->getGenotypeQuality();
						string zygosity = variantZygosity2str(variant_in_sample->getGenotype());

						ostringstream strs1;
						strs1 << zygosity << "\t";
						strs1 << GQ << "\t";

						if ( CRunSettings::getInstance()->get_use_genotyper_reads() ) {
							int genotyperReadsRef = variant_in_sample->getGenotyperReadsRef();
							int genotyperReadsVar = variant_in_sample->getGenotyperReadsVar();
							int genotyperReads = variant_in_sample->getGenotyperReads();
							double genotyperMutationFreq = variant_in_sample->getGenotyperMutationFreq();

							strs1 << genotyperReadsRef << "\t";
							strs1 << genotyperReadsVar << "\t";
							strs1 << genotyperReads << "\t";
							strs1 << genotyperMutationFreq << "\t";

						}
						else {
							int totalReadsRefMin = variant_in_sample->getTotalReadsRefMinus();
							int totalReadsRefPlus = variant_in_sample->getTotalReadsRefPlus();
							int totalReadsVarMin = variant_in_sample->getTotalReadsVarMinus();
							int totalReadsVarPlus = variant_in_sample->getTotalReadsVarPlus();
							int totalReadsVar = variant_in_sample->getTotalReadsVar();
							int totalReads = variant_in_sample->getTotalReads();
							double totalMutationFreq = variant_in_sample->getTotalMutationFreq();

							strs1 << totalReadsRefMin << "\t" ;
							strs1 << totalReadsRefPlus << "\t" ;
							strs1 << totalReadsVarMin << "\t" ;
							strs1 << totalReadsVarPlus << "\t" ;
							strs1 << totalReadsVar << "\t" ;
							strs1 << totalReads << "\t";
							strs1 << totalMutationFreq << "\t";
						}

// TODO: fetch extra format fields from variant in sample
//						const vector<string>* extra_format_fields = CRunSettings::getInstance()->get_extra_format_fields();


						const vector<string>* extra_format_fields = CRunSettings::getInstance()->get_extra_format_fields();
						for (vector<string>::const_iterator extra_format_fields_it = extra_format_fields->begin();
								extra_format_fields_it != extra_format_fields->end() ;
								extra_format_fields_it++ ) {
							if ( variant_in_sample->getFormat(*extra_format_fields_it).empty() ) {
								strs1 << "\t";
							}
							else {
								strs1 << variant_in_sample->getFormat(*extra_format_fields_it) << "\t";
							}
						}


//						for(map<string, string >::const_iterator it = variant_in_sample->getAllFormats()->begin(); it != variant_in_sample->getAllFormats()->end(); it++)
//						{
//						    string key = it->first;
//						    string value = it->second;
//						    //Do something
//						    //cout << "parsing " << key << ":" << value << endl;
//						    strs1 << value << "\t";
//						}


						family_members_info.append(strs1.str());
					}

				}
				if ( !found_individual ) {
					if (CRunSettings::getInstance()->get_use_genotyper_reads()) {
						family_members_info.append("\t").append("\t").append("\t").append("\t").append("\t").append("\t");
					}
					else {
						family_members_info.append("\t").append("\t").append("\t").append("\t").append("\t").append("\t").append("\t").append("\t").append("\t");
					}
					for ( unsigned int i=0; i<CRunSettings::getInstance()->get_extra_format_fields()->size(); i++) {
						family_members_info.append("\t");
					}
				}
			}

			/* Printing individual variant info */
			stream << (family_members_str.substr(0, family_members_str.size()-1)) << "\t";
			stream << family_members_info;
			stream << endl;

			printedFamilies.insert(sample->getFamily());
		}
	}
	return stream;

}

