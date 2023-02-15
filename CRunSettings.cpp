/*
 * CGlobal.cpp
 *
 *  Created on: 2011-10-31
 *  Author: Alexandre Dionne-Laporte
 */

#include "CRunSettings.h"
#include "ITypes.h"

#include <fstream>
#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace boost;
using namespace boost::program_options;
using namespace std;

/* helper function */
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
    return os;
}

const char* CRunSettings::VERSION = "Segregation Analysis v3.3.0 (supporting phased genotypes)";
const int CRunSettings::VERBOSE_SILENT = 0;
const int CRunSettings::VERBOSE_LOG = 1;
const int CRunSettings::VERBOSE_ERROR = 2;

const int CRunSettings::EXIT_OK_FINISH=0;
const int CRunSettings::EXIT_OK_HELP=0;
const int CRunSettings::EXIT_FAIL_initializeRun=1;
const int CRunSettings::EXIT_FAIL_VCF_INIT_ERROR=2;


/** Global static pointer used to ensure a single instance of the class. */
CRunSettings* CRunSettings::_instance = new CRunSettings();

void CRunSettings::initializeRun ( int ac, char *av[] ) throw (ERunSettings) {

	/*
	 * Declare an options description instance which will include all the options
	 */
	variables_map vm;
	vector<string> unrecognized_options;

	// Declare three groups of options.
	options_description general_options("General options", 200);
	general_options.add_options()
		( "help,h",
				value< bool >(&_help)->default_value(false)->implicit_value(true)->zero_tokens(),
				"produce an help message and exit" )
		( "verbose",
				value<int>(&_verbose_level)->default_value(0)->implicit_value(1),
				"ouput more details [0,1,2]" )
		( "version",
				value< bool >(&_version)->default_value(false)->implicit_value(true)->zero_tokens(),
				"print program version and details and exit" )
		;

	//Default values format filters
//	_format_filters->push_back("ADS");
//	_format_filters->push_back("AD");
//	_format_filters->push_back("DP");
//	string format_filters_help_str="";
//	for ( vector<string>::iterator i=_format_filters->begin(); i != _format_filters->end(); i++ ) {
//		format_filters_help_str.append(*i).append(" ");
//	}

	// variant filtration
	options_description seg_options("Treshold values for variant filtration", 200);
	seg_options.add_options()
		( "proband-var-reads",
				value< int >(&_min_proband_variant_reads)->default_value(1),
				"min. variant reads in probands"  )
		( "proband-coverage",
				value< int >(&_min_proband_coverage)->default_value(1),
				"min. coverage in probands"  )
		( "proband-mutfreq",
				value< double >(&_min_proband_mutfreq)->default_value(0),
				"min. mutation frequency in probands"  )
		( "proband-genotype-quality",
				value< int >(&_min_proband_genotype_quality)->default_value(0),
				"min. genotype quality in probands"  )
		( "familial-control-var-reads",
				value< int >(&_min_familial_control_variant_reads)->default_value(1),
				"min. variant reads in familial controls"  )
		( "familial-control-coverage",
				value< int >(&_min_familial_control_coverage)->default_value(1),
				"min. coverage in familial controls"  )
		( "familial-control-mutfreq",
				value< double >(&_min_familial_control_mutfreq)->default_value(0),
				"min. mutation frequency in familial controls"  )
		( "familial-control-genotype-quality",
				value< int >(&_min_familial_control_genotype_quality)->default_value(0),
				"min. genotype quality in familial controls"  )
		( "external-control-var-reads",
				value< int >(&_min_external_control_variant_reads)->default_value(3),
				"min. variant reads in external controls"  )
		( "external-control-coverage",
				value< int >(&_min_external_control_coverage)->default_value(6),
				"min. coverage in external controls"  )
		( "external-control-mutfreq",
				value< double >(&_min_external_control_mutfreq)->default_value(((double)(0.2)),"0.2"),
				"min. mutation frequency in external controls"  )
		( "external-control-genotype-quality",
				value< int >(&_min_external_control_genotype_quality)->default_value(20),
				"min. genotype quality in external controls"  )
		( "show-probands-with-variants",
				value< int >(&_show_probands_with_variants)->default_value(10),
				"max. proband names to show" )
		( "show-familial-controls-with-variants",
				value< int >(&_show_familial_controls_with_variants)->default_value(10),
				"max. familial control names to show" )
		( "show-external-controls-with-variants",
				value< int >(&_show_external_controls_with_variants)->default_value(10),
				"max. external control names to show" )
		( "use-genotyper-reads",
				value< bool >(&_use_genotyper_reads)->default_value(false)->implicit_value(true)->zero_tokens(),
				"use AD reads (instead of ADS) for coverage filtering" )
//		( "filter-format",
//				value< vector<string> >(_format_filters)->default_value(*_format_filters, format_filters_help_str)->multitoken(),
//				"which field in format should we use for filtration")
//		( "lenient-parsing",
//				value< bool >(&_lenient_parsing)->default_value(false)->implicit_value(true)->zero_tokens(),
//				"very relaxed VCF parsing" )
		( "use-indel-filtering",
				value< bool >(&_use_indel_filtering)->default_value(false)->implicit_value(true)->zero_tokens(),
				"use filters for indel variants" )
		( "use-format-low-stringency",
				value< bool >(&_format_low_stringency)->default_value(false)->implicit_value(true)->zero_tokens(),
				"use a low format stringency for parsing (implies use-genotyper-reads)" )
		( "show-all-variants",
				value< bool >(&_show_all_variants)->default_value(false)->implicit_value(true)->zero_tokens(),
				"show wildtype calls in probands and variants found in controls" )
//		( "show-all-variants-in-family",
//				value< bool >(&_show_all_variants)->default_value(false)->implicit_value(true)->zero_tokens(),
//				"show wildtype calls in probands and variants found in familial controls" )
		;


	//Default values for the info fields
//	_extra_info_fields->push_back("VC");
//	_extra_info_fields->push_back("VFT");
//	_extra_info_fields->push_back("VT");
//	_extra_info_fields->push_back("Gene");
//	_extra_info_fields->push_back("1KG");
//	_extra_info_fields->push_back("CG");
//	_extra_info_fields->push_back("DA");
//	_extra_info_fields->push_back("ESP6500");
//	_extra_info_fields->push_back("GERP++");
//	_extra_info_fields->push_back("LRT");
//	_extra_info_fields->push_back("MT");
//	_extra_info_fields->push_back("NCdbSNP_rs");
//	_extra_info_fields->push_back("PhyloP");
//	_extra_info_fields->push_back("PPv2");
//	_extra_info_fields->push_back("SIFT");

	string info_fields_help_str="";
	for ( vector<string>::iterator i=_extra_info_fields->begin(); i != _extra_info_fields->end(); i++ ) {
		info_fields_help_str.append(*i).append(" ");
	}

	//Default values for the format fields
//	_extra_format_fields->push_back("ADS");
//	_extra_format_fields->push_back("AD");
//	_extra_format_fields->push_back("DP");
	string format_fields_help_str="";
	for ( vector<string>::iterator i=_extra_format_fields->begin(); i != _extra_format_fields->end(); i++ ) {
		format_fields_help_str.append(*i).append(" ");
	}

	options_description file_options("Input/output parameters", 200);
	file_options.add_options()
		("vcf-file,v", value< vector<string> >(&_vcf_input_files)->multitoken(), "vcf input file(s) and optional filters. \nFormat: file.vcf[:pR,pX,p%,pQ:fcR,fcX,fc%,fcQ:nfcR,nfcX,nfc%,ncfQ]")
		("ped-file,p", value<string>(&_pedigree_input_file), "pedigree input file")
		("out-file,o", value<string>(&_output_file), "output file")
		("add-info,i", value< vector<string> >(_extra_info_fields)->default_value(*_extra_info_fields, info_fields_help_str)->multitoken(), "supplementary INFO field(s) to output")
		("add-format,f", value< vector<string> >(_extra_format_fields)->default_value(*_extra_format_fields, format_fields_help_str)->multitoken(), "supplementary FORMAT field(s) to output")
		;

	// Add options to the common container
	_all->add(file_options).add(seg_options).add(general_options);
	parsed_options parsed = command_line_parser(ac, av).options(*_all).allow_unregistered().run();
	store(parsed, vm);
	notify(vm);
	unrecognized_options = collect_unrecognized(parsed.options, include_positional);

	/*
	 * General options
	 */
	if ( !_help && !_version ) {

		/* Verifying unrecognized options */
		if ( unrecognized_options.size() > 0 ) {
			string error("Unknown option(s): ");
			for (unsigned int i=0; i < unrecognized_options.size(); i++) {
				error.append("\"").append(unrecognized_options[i]).append("\" ");
			}
			throw ERunSettings(error, __FILE__, __LINE__);
		}

		/* VCF files option */
		if (!_vcf_input_files.empty())
		{
			for ( vector<string>::iterator i=_vcf_input_files.begin(); i != _vcf_input_files.end(); i++ ) {
				string vfc_file_arg = *i;
				string vfc_file = vfc_file_arg.substr(0, vfc_file_arg.find(':') );
				ifstream ifile( vfc_file.c_str() );
				if (!ifile) {
					throw ERunSettings("Input error: must provide a valid vcf file\n" + *i, __FILE__, __LINE__);
				}
			}
		}
		else {
			throw ERunSettings("Input error: must provide a vcf file", __FILE__, __LINE__);
		}

		/* PED file option */
		if ( !_pedigree_input_file.empty() )
		{
			ifstream ifile (_pedigree_input_file.c_str());
			if (!ifile) {
				throw ERunSettings("Input error: must provide a valid ped file\n" + _pedigree_input_file, __FILE__, __LINE__);
			}
		}
		else {
			throw ERunSettings("Input error: must provide a ped file\n" + _pedigree_input_file, __FILE__, __LINE__);
		}

		/*
		 * OUPUT file option
		 */
		if ( !_output_file.empty() )
		{
			ifstream ifile(_output_file.c_str());
		}
		else {
			throw ERunSettings("Input error: must provide an output file", __FILE__, __LINE__);
		}

		/*
		 * Special case with format_low_stringency
		 */
		if ( _format_low_stringency ) {
			_use_genotyper_reads = true;
		}

		/*
		 * Special case with use_genotyper_reads
		 */
//		_extra_format_fields->insert(_extra_format_fields->begin(), "DP");
//		_extra_format_fields->insert(_extra_format_fields->begin(), "AD");
//		if ( ! _use_genotyper_reads ) {
//			_extra_format_fields->insert(_extra_format_fields->begin(), "ADS");
//		}
		/*
		 * Add the default FORMATs at the beginning
		 */

	}

	/* Setting all required INFO FIELDS from the vcf */
	{
//		_required_info_fields->push_back("1KG");
//		_required_info_fields->push_back("CG");
//		_required_info_fields->push_back("DA");
//		_required_info_fields->push_back("ESP6500");
//		_required_info_fields->push_back("Gene");
//		_required_info_fields->push_back("GERP++");
//		_required_info_fields->push_back("LRT");
//		_required_info_fields->push_back("MT");
//		_required_info_fields->push_back("NCdbSNP_rs");
//		_required_info_fields->push_back("PhyloP");
//		_required_info_fields->push_back("PPv2");
//		_required_info_fields->push_back("SIFT");
//		_required_info_fields->push_back("VC");
//		_required_info_fields->push_back("VFT");
//		_required_info_fields->push_back("VT");
		for ( vector<string>::iterator it=_extra_info_fields->begin(); it != _extra_info_fields->end(); it++ ) {
			string val = *it;
			_required_info_fields->push_back(val);
		}
	}

	/* Setting all required FORMAT FIELDS from the vcf */
	{
		for ( vector<string>::iterator it=_extra_format_fields->begin(); it != _extra_format_fields->end(); it++ ) {
			string val = *it;
			_required_format_fields->push_back(val);
		}
	}

	_isInitialized = true;
}

string CRunSettings::getRunArguments() {

	std::stringstream ss;

	ss << output_seperator
	   << endl << CRunSettings::VERSION << endl << endl;

	ss << output_seperator << endl;
	ss << "RUN SETTINGS" << endl;
	for ( unsigned int i=0; i<_vcf_input_files.size(); i++ ) {
		ss << "vcf input file (" << i+1 << ")               \"" << _vcf_input_files[i] << "\"" << endl;
	}
	ss << "ped input file                   \"" << _pedigree_input_file << "\"" << endl;
	ss << "output file                      \"" << _output_file << "\"" << endl;
	for ( unsigned int i=0; i<_extra_info_fields->size(); i++ ) {
		ss << "supplementary info field (" << i+1 << ")     \"" << (*_extra_info_fields)[i] << "\"" << endl;
	}
	for ( unsigned int i=0; i<_extra_format_fields->size(); i++ ) {
		ss << "supplementary format field (" << i+1 << ")   \"" << (*_extra_format_fields)[i] << "\"" << endl;
	}
	ss << "proband variant reads               " << _min_proband_variant_reads << endl;
	ss << "fam. control variant reads          " << _min_familial_control_variant_reads << endl;
	ss << "external control variant reads      " << _min_external_control_variant_reads << endl;
	ss << "proband coverage                    " << _min_proband_coverage << endl;
	ss << "fam. control coverage               " << _min_familial_control_coverage << endl;
	ss << "external control coverage           " << _min_external_control_coverage << endl;
	ss << "proband mutation frequency          " << _min_proband_mutfreq << endl;
	ss << "fam. control mutation frequency     " << _min_familial_control_mutfreq << endl;
	ss << "external control mutation frequency " << _min_external_control_mutfreq << endl;
	ss << "proband genotype quality            " << _min_proband_genotype_quality << endl;
	ss << "fam. control genotype quality       " << _min_familial_control_genotype_quality << endl;
	ss << "external control genotype quality   " << _min_external_control_genotype_quality << endl;
	ss << "use genotyper reads counts filters  " << _use_genotyper_reads << endl;
//	ss << "use relaxed VCF parsing             " << _lenient_parsing << endl;
	ss << "use filters in indels               " << _use_indel_filtering << endl;
	ss << "use low format stringency           " << _format_low_stringency << endl;
	ss << "show external controls w/ variant   " << _show_external_controls_with_variants << endl;
	ss << "show filtered and control variants  " << _show_all_variants << endl;
	ss << "verbosity level                     " << _verbose_level << endl;

	return ss.str();
}

bool CRunSettings::get_is_required_info_fields(string field) {
	bool returnValue=false;
	for (vector<string>::iterator it=_required_info_fields->begin(); it!=_required_info_fields->end(); it++) {
		if ( it->compare(field.c_str()) == 0 ) {
			returnValue=true;
			break;
		}
	}
	return returnValue;
}

bool CRunSettings::get_is_required_format_fields(string field) {
	bool returnValue=false;
	for (vector<string>::iterator it=_required_format_fields->begin(); it!=_required_format_fields->end(); it++) {
		if ( it->compare(field.c_str()) == 0 ) {
			returnValue=true;
			break;
		}
	}
	return returnValue;
}

bool CRunSettings::insert_info_description (string name, string description) {
	bool returnValue=false;
	if ( _allInfoDescriptions->find(name) == _allInfoDescriptions->end() ) {
		(*_allInfoDescriptions)[name] = description;
		returnValue=true;
	}
	return returnValue;
}

bool CRunSettings::insert_format_description (string name, string description) {
	bool returnValue=false;
	if ( _allFormatDescriptions->find(name) == _allFormatDescriptions->end() ) {
		(*_allFormatDescriptions)[name] = description;
		returnValue=true;
	}
	return returnValue;
}

CRunSettings::~CRunSettings() {
	/* VERBOSE output */
	if (CRunSettings::getInstance()->get_verbose_level() > CRunSettings::VERBOSE_LOG)
		cerr << "CVariant::~CVariant" << endl;
	_extra_info_fields->clear();
	delete _extra_info_fields;
	_required_info_fields->clear();
	delete _required_info_fields;
	_extra_format_fields->clear();
	delete _extra_format_fields;
	_required_format_fields->clear();
	delete _required_format_fields;
	_allInfoDescriptions->clear();
	delete _allInfoDescriptions;
};
