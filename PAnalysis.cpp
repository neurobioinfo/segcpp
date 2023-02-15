/*
 * PAnalysis.cpp
 *
 *  Created on: 2011-11-09
 *  Author: Alexandre Dionne-Laporte
 */

#include "CRunSettings.h"
#include "CVariant.h"
#include "CSample.h"
#include "UPedFileParser.h"
#include "CVCFFile.h"
#include "CChromosome.h"
#include "CFamily.h"
#include "ITypes.h"
#include "UStringUtils.h"

#include <iostream>
#include <fstream>
#include <iterator>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

using namespace boost;
using namespace boost::program_options;

/**
 * Helper function
 */
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
	copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
	return os;
}

void printFIELDS(CVCFFile& vcfFile) {
	cout << endl << output_seperator << endl;
	cout << "VCF INFO FIELDS DESCRIPTION" << endl;
	map<string, VariantFieldType> infoTypes = vcfFile.getInfoTypes();
	map<string, int> infoCounts = vcfFile.getInfoCounts();
	map<string, string> infoDescriptions = vcfFile.getInfoDescriptions();
	for (std::map<string, string>::const_iterator it = infoDescriptions.begin(); it
	!= infoDescriptions.end(); it++) {
		cout << it->first << "\t" << it->second << "\t"
				<< variantFieldNumber2str(infoCounts[it->first]) << "\t"
				<< variantFieldType2str(infoTypes[it->first]) << "\t" << endl;
	}
	cout << endl << output_seperator << endl;
	cout << "VCF FORMAT FIELDS DESCRIPTIONS" << endl;
	map<string, VariantFieldType> formatTypes = vcfFile.getFormatTypes();
	map<string, int> formatCounts = vcfFile.getFormatCounts();
	map<string, string> formatDescriptions = vcfFile.getFormatDescriptions();
	for (std::map<string, string>::const_iterator it =
			formatDescriptions.begin(); it != formatDescriptions.end(); it++) {
		cout << it->first << "\t" << it->second << "\t"
				<< variantFieldNumber2str(formatCounts[it->first]) << "\t"
				<< variantFieldType2str(formatTypes[it->first]) << "\t" << endl;
	}
	cout << endl << output_seperator << endl;
	cout << "VCF SAMPLES " << endl << vcfFile.getSampleNames() << endl;
}

string getHeader(vector<string> extraInfoFieldDescriptions, vector<string> extraFormatFieldDescriptions) {

	string header = string("Family ID\t"
			"Chr\tPosition\t"
			"Reference Allele\t"
			"Mutant Allele\t");

	header += /*CRunSettings::getInstance()->get_info_field_description(CVCFFile::VARIANT_CLASS) + "\t"
			+ CRunSettings::getInstance()->get_info_field_description(CVCFFile::GENE) + "\t"
			+ CRunSettings::getInstance()->get_info_field_description(CVCFFile::VARIANT_FUNCTION_TYPE) + "\t"
			+ CRunSettings::getInstance()->get_info_field_description(CVCFFile::VARIANT_TYPE) + "\t"
			+ CRunSettings::getInstance()->get_info_field_description(CVCFFile::DETAILED_ANNOTATION) + "\t"
			+ */"Number of probands w/ variant in gene\t";// +
	//			CRunSettings::getInstance()->get_info_field_description(CVCFFile::C_DBSNP_RS) + "\t" +
	//			CRunSettings::getInstance()->get_info_field_description(CVCFFile::NC_DBSNP_RS) + "\t" +
	//			CRunSettings::getInstance()->get_info_field_description(CVCFFile::THOUSAND) + "\t" +
	//			CRunSettings::getInstance()->get_info_field_description(CVCFFile::CG) + "\t" +
	//			CRunSettings::getInstance()->get_info_field_description(CVCFFile::EVS) + "\t" +
	//			CRunSettings::getInstance()->get_info_field_description(CVCFFile::SIFT) + "\t" +
	//			CRunSettings::getInstance()->get_info_field_description(CVCFFile::PPV2) + "\t" +
	//			CRunSettings::getInstance()->get_info_field_description(CVCFFile::LRT) + "\t" +
	//			CRunSettings::getInstance()->get_info_field_description(CVCFFile::MT) + "\t" +
	//			CRunSettings::getInstance()->get_info_field_description(CVCFFile::PHYLOP) + "\t" +
	//			CRunSettings::getInstance()->get_info_field_description(CVCFFile::GERP) + "\t";

	for (vector<string>::iterator i = extraInfoFieldDescriptions.begin(); i
	!= extraInfoFieldDescriptions.end(); i++) {
		header += (*i).append("\t");
	}

	header += //"Freeze Set Allele Frequency\t"
			//			"Dels\t"
			//			"Strand Bias\t"
			//			"Mapping Quality\t"
			//			"Quality Depth\t"
			"Freeze Set Filter\t";
	//			"Cultprit\t";

	header += "Segregation Allele Count\t"
			"Segregation Allele Frequency\t";

	//	if ( CRunSettings::getInstance()->get_show_external_controls_with_variants() > 0 ) {
	header += "Probands w/ variant individual names\t";
	header += "Familial controls w/ variant individual names\t";
	header += "External controls w/ variant individual names\t";
	//	}


	header += "Total subjects filtered\t"
			"Total subjects wildtype\t"
			"Total subjects w/ variant\t"
			"Total affected w/ variant\t"
			"Total controls w/ variant\t"
			"External controls filtered\t"
			"External controls wildtype\t"
			"External controls w/ variant\t"
			"External controls w/ homozygote variant\t"
			"Family subjects w/ variant\t"
			"Family affected filtered\t"
			"Family affected wildtype\t"
			"Family affected w/ variant\t"
			"Family controls filtered\t"
			"Family controls wildtype\t"
			"Family controls w/ variant\t";

	if (CRunSettings::getInstance()->get_use_genotyper_reads()) {
		header += "Family members\t"
				"Family member zygosity\t"
				"Family member genotype qual\t"
				//				"Family member ref- cov\t"
				//				"Family member ref+ cov\t"
				//				"Family member var- cov\t"
				//				"Family member var+ cov\t"
				//				"Family member var cov\t"
				//				"Family member total cov\t"
				//				"Family member reads mut freq\t";
				"Family member genotyper ref cov\t"
				"Family member genotyper var cov\t"
				"Family member genotyper cov\t"
				"Family member genotyper mut freq\t";
	} else {
		header += "Family members\t"
				"Family member zygosity\t"
				"Family member genotype qual\t"
				"Family member ref- cov\t"
				"Family member ref+ cov\t"
				"Family member var- cov\t"
				"Family member var+ cov\t"
				"Family member var cov\t"
				"Family member total cov\t"
				"Family member reads mut freq\t";
		//			"Family member genotyper ref cov\t"
		//			"Family member genotyper var cov\t"
		//			"Family member genotyper cov\t"
		//			"Family member genotyper mut freq";
	}

	for (vector<string>::iterator i = extraFormatFieldDescriptions.begin(); i
	!= extraFormatFieldDescriptions.end(); i++) {
		header += (*i).append("\t");
	}

	return header;
}

int main(int ac, char *av[]) {
	CRunSettings* runSettings = CRunSettings::getInstance();

	/* Parsing run's arguments */
	try {
		runSettings->initializeRun(ac, av);
	} catch (boost::program_options::error& e) {
		cerr << output_seperator << endl << e.what() << endl
				<< output_seperator << endl;
		return CRunSettings::EXIT_FAIL_initializeRun;
	} catch (ERunSettings& e) {
		cerr << e.message() << endl;
		if (runSettings->get_options_description() != NULL) {
			cerr << output_seperator << endl
					<< *(runSettings->get_options_description());
		}
		cerr << output_seperator << endl;
		return CRunSettings::EXIT_FAIL_initializeRun;
	}
	/* HELP option */
	if (runSettings->get_help()) {
		cout << output_seperator << endl
			 << (*(runSettings->get_options_description()))
			 << output_seperator << endl;
		return CRunSettings::EXIT_OK_HELP;
	}
	/* VERSION option */
	if (runSettings->get_version()) {
		cout << output_seperator << endl << CRunSettings::VERSION << endl
			 << output_seperator << endl;
		return CRunSettings::EXIT_OK_HELP;
	}
	/* REGULAR output */
	cout << runSettings->getRunArguments();

	/* Parsing Ped file */
	ifstream ped_stream(string(runSettings->get_pedigree_input_file()).c_str());
	vector<CSample*> * samples = UPedFileParser::parseFile(ped_stream);

	/* Parsing VCF files */
	vector<CVariant*>* variants = new vector<CVariant*> ();
	int total_variants = 0;
	vector<string> vcf_input_files = runSettings->get_vcf_input_files();
	const vector<string>* extra_info_fields =runSettings->get_extra_info_fields();
	vector<string> extra_info_fields_description;
	const vector<string>* extra_format_fields =runSettings->get_extra_format_fields();
	vector<string> extra_format_fields_description;

	/* Default filters */
	int min_proband_variant_reads =	runSettings->get_min_proband_variant_reads();
	int min_familial_control_variant_reads =runSettings->get_min_familial_control_variant_reads();
	int min_external_control_variant_reads =runSettings->get_min_external_control_variant_reads();
	int min_proband_coverage = runSettings->get_min_proband_coverage();
	int min_familial_control_coverage =	runSettings->get_min_familial_control_coverage();
	int min_external_control_coverage =	runSettings->get_min_external_control_coverage();
	double min_proband_mutfreq = runSettings->get_min_proband_mutfreq();
	double min_familial_control_mutfreq =runSettings->get_min_familial_control_mutfreq();
	double min_external_control_mutfreq =runSettings->get_min_external_control_mutfreq();
	int min_proband_genotype_quality =runSettings->get_min_proband_genotype_quality();
	int min_familial_control_genotype_quality =runSettings->get_min_familial_control_genotype_quality();
	int min_external_control_genotype_quality =runSettings->get_min_external_control_genotype_quality();

	for (vector<string>::iterator i = vcf_input_files.begin(); i != vcf_input_files.end(); i++) {
		CVCFFile *vcfFile = new CVCFFile();
		string file_name = *i;
		string vcf_input_text = *i;

		/* Specific filters for each input file */
		char_separator<char> vcf_input_separator(":");
		tokenizer<char_separator<char> > vcf_input_tokens(vcf_input_text,
				vcf_input_separator);
		int vcf_input_num = 0;
		BOOST_FOREACH(string status_filter_text, vcf_input_tokens)
		{
			vcf_input_num++;
			char_separator<char> status_filter_separator(",");
			tokenizer<char_separator<char> > status_filter_token(
					status_filter_text, status_filter_separator);
			int status_filter_num = 0;
			BOOST_FOREACH(string status_filter_sub_text, status_filter_token)
			{
				status_filter_num++;

				switch (vcf_input_num) {

				case 1:
					// This is the vcf file itself
					file_name = status_filter_text;
					break;
				case 2:
					// Process the proband options
					switch (status_filter_num) {
					case 1:
						// Process the prob_var_reads
						min_proband_variant_reads= atoi(status_filter_sub_text.c_str());
						break;
					case 2:
						// Process the prob_cov
						min_proband_coverage
						= atoi(status_filter_sub_text.c_str());
						break;
					case 3:
						// Process the prob_mutfreq
						min_proband_mutfreq
						= atof(status_filter_sub_text.c_str());
						break;
					case 4:
						// Process the prob genotype_quality
						min_proband_genotype_quality
						= atoi(status_filter_sub_text.c_str());
						break;
					default:
						// Process for all other cases.
						cerr << "No such option, see "
						<< __FILE__
						<< " at line "
						<< __LINE__ << endl;
					}
					break;

					case 3:
						// Process the familial controls options
						switch (status_filter_num) {
						case 1:
							// Process the fam_cont_var_reads
							min_familial_control_variant_reads
							= atoi(status_filter_sub_text.c_str());
							break;
						case 2:
							// Process the fam_cont_cov
							min_familial_control_coverage
							= atoi(status_filter_sub_text.c_str());
							break;
						case 3:
							// Process the fam_cont_mutfreq
							min_familial_control_mutfreq
							= atof(status_filter_sub_text.c_str());
							break;
						case 4:
							// Process the fam_cont_genotype_quality
							min_familial_control_genotype_quality
							= atoi(status_filter_sub_text.c_str());
							break;
						default:
							// Process for all other cases.
							cerr << "No such option, see "
							<< __FILE__
							<< " at line "
							<< __LINE__ << endl;
						}
						break;

						case 4:
							// Process the external controls options
							switch (status_filter_num) {
							case 1:
								// Process the non_fam_cont_var_reads
								min_external_control_variant_reads
								= atoi(status_filter_sub_text.c_str());
								break;
							case 2:
								// Process the non_fam_cont_cov
								min_external_control_coverage
								= atoi(status_filter_sub_text.c_str());
								break;
							case 3:
								// Process the non_fam_cont_mutfreq
								min_external_control_mutfreq
								= atof(status_filter_sub_text.c_str());
								break;
							case 4:
								// Process the non_fam_cont_genotype_quality
								min_external_control_genotype_quality
								= atoi(status_filter_sub_text.c_str());
								break;
							default:
								// Process for all other cases.
								cerr << "No such option, see "
								<< __FILE__
								<< " at line "
								<< __LINE__ << endl;
							}
							break;

							default:
								// Process for all other cases.
								cerr << "No such option, see "
								<< __FILE__ << " at line "
								<< __LINE__ << endl;
				}
			}
		}

		/* opening vcf file and check header conformity */
		try {
			vcfFile->openFile(file_name);

			/* CHECK if sample defined in another vcf file */
			vector<string> sample_names = vcfFile->getSampleNames();
			for (vector<string>::iterator it = sample_names.begin(); it
			!= sample_names.end(); it++) {
				CSample* current_sample = CSample::getSample(*it);
				if (current_sample != NULL && current_sample->getVariantsInSample()->size() > 0) {
					cerr << "Duplicate sample: already have "
						 << current_sample->getVariantsInSample()->size()
						 << " variants for sample "
						 << current_sample->getName() << endl;
				}
			}

			/* Verifying that each extra INFO field is present */
			{
				for (vector<string>::const_iterator i = extra_info_fields->begin(); i != extra_info_fields->end(); i++) {
					map<string, string> infoDescriptions = vcfFile->getInfoDescriptions();
					if (infoDescriptions.find(*i) == infoDescriptions.end()) {
						cerr << endl << output_seperator << endl;
						cerr << "Invalid InfoField for VCF: \"" << *i << "\""
								<< endl;
						return CRunSettings::EXIT_FAIL_VCF_INIT_ERROR;
					} else {
						bool foundDesc = false;
						for (vector<string>::iterator eifd = extra_info_fields_description.begin(); eifd != extra_info_fields_description.end(); eifd++) {
							string infoDescCompare = infoDescriptions[*i];
							infoDescCompare = infoDescCompare.append(" (").append(*i).append(")");
							if ((*eifd).compare(infoDescCompare) == 0) { foundDesc = true;	}
//							if ( infoDescriptions.find(*i) != infoDescriptions.end() ) { foundDesc = true;	}
						}
						if (!foundDesc) {
							string toInsert = infoDescriptions[*i];
							toInsert = toInsert.append(" (").append(*i).append(")");
							extra_info_fields_description.push_back(toInsert);
						}
					}
				}
				/** Add all info fields description in a map */
				map<string, string> infoDescriptions = vcfFile->getInfoDescriptions();
				for (map<string, string>::iterator iter = infoDescriptions.begin(); iter != infoDescriptions.end(); iter++) {
					CRunSettings::getInstance()->insert_info_description(iter->first, iter->second);
				}
			}

			/* Verifying that each extra FORMAT field is present */
			{
				vector<string> needed_format_fields;
//We need AD and DP to output sample specific values
//				needed_format_fields.push_back("GQ");
				needed_format_fields.push_back("DP");
				needed_format_fields.push_back("AD");
				if ( !CRunSettings::getInstance()->get_format_low_stringency() ) {
					needed_format_fields.push_back("ADS");
				}
				for (vector<string>::const_iterator i = needed_format_fields.begin(); i != needed_format_fields.end(); i++) {
					map<string, string> formatDescriptions = vcfFile->getFormatDescriptions();
					if (formatDescriptions.find(*i) == formatDescriptions.end()) {
						cerr << endl << output_seperator << endl;
						cerr << "FormatField is not in VCF: \"" << *i << "\""
							 << endl;
						return CRunSettings::EXIT_FAIL_VCF_INIT_ERROR;
					}
				}
				for (vector<string>::const_iterator i = extra_format_fields->begin(); i != extra_format_fields->end(); i++) {
					map<string, string> formatDescriptions = vcfFile->getFormatDescriptions();
					if (formatDescriptions.find(*i) == formatDescriptions.end()) {
						cerr << endl << output_seperator << endl;
						cerr << "FormatField is not in VCF: \"" << *i << "\""
							 << endl;
						return CRunSettings::EXIT_FAIL_VCF_INIT_ERROR;
					} else {
						bool foundDesc = false;
						for (vector<string>::iterator effd = extra_format_fields_description.begin(); effd != extra_format_fields_description.end(); effd++) {
							string formatDescCompare = formatDescriptions[*i];
							formatDescCompare = formatDescCompare.append(" (").append(*i).append(")");
							if ((*effd).compare(formatDescCompare) == 0) { foundDesc = true;	}
						}
						if (!foundDesc) {
							string toInsert = formatDescriptions[*i];
							toInsert = toInsert.append(" (").append(*i).append(")");
							extra_format_fields_description.push_back(toInsert);
						}
					}
				}

				/** Add all format fields description in a map */
				map<string, string> formatDescriptions = vcfFile->getFormatDescriptions();
				for (map<string, string>::iterator iter = formatDescriptions.begin(); iter != formatDescriptions.end(); iter++) {
					CRunSettings::getInstance()->insert_format_description(iter->first, iter->second);
				}
			}

			/* HIGH VERBOSE output */
			if (runSettings->get_verbose_level() > CRunSettings::VERBOSE_LOG)
				printFIELDS(*vcfFile);
		} catch (EVCFFile& e) {
			cerr << endl << output_seperator << endl << e.message() << endl << output_seperator << endl;
//			printFIELDS(*vcfFile);
			return CRunSettings::EXIT_FAIL_VCF_INIT_ERROR;
		}

		/* Parsing VCF file content */
		vector<CVariant*>* vcf_variants = vcfFile->parseFile(
				min_proband_variant_reads, min_proband_coverage,
				min_proband_mutfreq, min_proband_genotype_quality,
				min_familial_control_variant_reads,
				min_familial_control_coverage, min_familial_control_mutfreq,
				min_familial_control_genotype_quality,
				min_external_control_variant_reads,
				min_external_control_coverage, min_external_control_mutfreq,
				min_external_control_genotype_quality);

		//Adding to total
		variants->insert(variants->begin(), vcf_variants->begin(), vcf_variants->end());
		total_variants += (vcfFile->getLineNumber());

		delete vcfFile;
	}

	/* Verifying that we got all samples from the PED file */
	int samples_without_vairants = 0;
	vector<CSample*> all_samples = CSample::getSamples();
	for (vector<CSample*>::iterator it = all_samples.begin(); it != all_samples.end(); ++it) {
		CSample* current_sample = *it;
		if (current_sample->getVariantsInSample()->empty()) {
			samples_without_vairants++;
		}
	}

	/* VERBOSE output */
	if (runSettings->get_verbose_level() > CRunSettings::VERBOSE_SILENT) {
		cout << endl << output_seperator << endl;
		cout << "PEDIGREE AND VARIANT INFO" << endl;
		for (vector<CSample*>::iterator i = samples->begin(); i	!= samples->end(); i++) {
			cout << **i << ", variants: " << (*i)->getVariantsInSample()->size() << endl;
		}
	}

	/* REGULAR output */
	cout << endl << output_seperator << endl;
	cout << "RUN STATS" << endl;
	cout << "number of families:            \t"
		 << CFamily::getAllFamilies()->size() << endl;
	cout << "number of samples defined:     \t" << samples->size() << endl;
	cout << "number of samples w/o variants:\t" << samples_without_vairants
		 << endl;
	cout << "total variant sites:           \t" << total_variants << endl;

	/* Generating output file */
	/* VERBOSE output */
	if (runSettings->get_verbose_level() > CRunSettings::VERBOSE_SILENT) {
		cout << endl << output_seperator << endl;
		cout << "GENERATING OUTPUT..." << flush;
	}

	/* File output here */
	ofstream output_stream(string(runSettings->get_output_file()).c_str());
	output_stream << getHeader(extra_info_fields_description, extra_format_fields_description) << endl;
	vector<CChromosome*> * chromosomes = CChromosome::getAllChromosomes();
	//loop over the chromosomes
	for (vector<CChromosome*>::iterator chr_it = chromosomes->begin(); chr_it < chromosomes->end(); chr_it++) {
		CChromosome * chromosome = *chr_it;
		const vector<CVariant*> * variants = chromosome->getAllVariants();
		//loop over each variants in the current chromosome
		for (vector<CVariant*>::const_iterator var_it = variants->begin(); var_it	< variants->end(); var_it++) {
			CVariant *variant = *var_it;
			output_stream << *variant;
		}
	}
	/* VERBOSE output */
	if (runSettings->get_verbose_level() > CRunSettings::VERBOSE_SILENT) {
		cout << "\tdone!" << endl;
	}

	return CRunSettings::EXIT_OK_FINISH;
}

