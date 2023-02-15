/*
 * CRunSettings.h
 *
 * Created on: 2011-10-31
 * Author: Alexandre Dionne-Laporte
 */

#ifndef CRUN_SETTINGS_H_
#define CRUN_SETTINGS_H_

#include <boost/program_options/options_description.hpp>
#include <string>
#include <iostream>

using namespace std;
using namespace boost::program_options;

class ERunSettings : public std::exception {
protected:
	string _msg;
	string _file;
	string _line;
public:
	ERunSettings(const string msg, const string file, const int line) : _msg(msg), _file(file) {
		std::stringstream ss;
		ss << line;
		_line = ss.str();
	};
	virtual ~ERunSettings() throw() {}
	const string message() {
		std::stringstream ss;
		ss << _msg << endl << "In file " << _file << ", line " << _line << endl;;
		return ss.str();
	}
//	virtual const char* what() const throw() {
//		string ret;
//		ret.append(_msg).append("\nIn file ").append(_file).append(", line ").append(_line);
//		cout << ret <<endl;
//		return ret.c_str();
//	}
};

class CRunSettings {

public:
	/**
	 * CONSTANTS
	 */
	static const char* VERSION;
	static const int VERBOSE_SILENT;
	static const int VERBOSE_LOG;
	static const int VERBOSE_ERROR;

	static const int EXIT_OK_FINISH;
	static const int EXIT_OK_HELP;
	static const int EXIT_FAIL_initializeRun;
	static const int EXIT_FAIL_VCF_INIT_ERROR;


	/**
	 * GETTERS / SETTERS
	 */
	inline const bool get_help() { return _help; }
	inline const bool get_version() { return _version; }
	inline const int get_verbose_level() { return _verbose_level; }

	inline const int get_min_proband_variant_reads() { return _min_proband_variant_reads; }
	inline const int get_min_familial_control_variant_reads() { return _min_familial_control_variant_reads; }
	inline const int get_min_external_control_variant_reads() { return _min_external_control_variant_reads; }
	inline const int get_min_proband_coverage() { return _min_proband_coverage; }
	inline const int get_min_familial_control_coverage() { return _min_familial_control_coverage; }
	inline const int get_min_external_control_coverage() { return _min_external_control_coverage; }
	inline const double get_min_proband_mutfreq() { return _min_proband_mutfreq; }
	inline const double get_min_familial_control_mutfreq() { return _min_familial_control_mutfreq; }
	inline const double get_min_external_control_mutfreq() { return _min_external_control_mutfreq; }
	inline const int get_min_proband_genotype_quality() { return _min_proband_genotype_quality; }
	inline const int get_min_familial_control_genotype_quality() { return _min_familial_control_genotype_quality; }
	inline const int get_min_external_control_genotype_quality() { return _min_external_control_genotype_quality; }

	inline const bool get_indel_filtering() { return _use_indel_filtering; }
	inline const bool get_format_low_stringency() { return _format_low_stringency; }
	inline const unsigned int get_show_probands_with_variants() { return _show_probands_with_variants; }
	inline const unsigned int get_show_familial_controls_with_variants() { return _show_familial_controls_with_variants; }
	inline const unsigned int get_show_external_controls_with_variants() { return _show_external_controls_with_variants; }
	inline const bool get_show_all_variants() { return _show_all_variants; }
	inline const bool get_use_genotyper_reads() { return _use_genotyper_reads; }
//	inline const bool get_lenient_parsing() { return _lenient_parsing; }

	inline const vector<string> get_vcf_input_files() { return _vcf_input_files; }
	inline const string get_pedigree_input_file() { return _pedigree_input_file; }
	inline const string get_output_file() { return _output_file; }
	inline const vector<string>* get_extra_info_fields() { return _extra_info_fields; }
	inline const vector<string>* get_extra_format_fields() { return _extra_format_fields; }
	inline const string get_info_field_description(string name) const { return (*_allInfoDescriptions)[name]; }

	inline const options_description* get_options_description() { return _all; }
	inline const  bool getIsInitialized() { return _isInitialized; }

	static CRunSettings* getInstance() { return _instance; }

	/**
	 * PUBLIC METHODS
	 */
	void initializeRun ( int ac, char *av[] ) throw (ERunSettings) ;
	string getRunArguments();
	bool get_is_required_info_fields(string field);
	bool get_is_required_format_fields(string field);
	bool insert_info_description (string name, string description);
	bool insert_format_description (string name, string description);

private:
	/**
	 * RUN VARIABLES
	 */

	/* general options */
	bool _help;
	bool _version;
	int _verbose_level;

	/* filtration options */
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
	bool _use_indel_filtering;
	bool _format_low_stringency;
	int _show_probands_with_variants;
	int _show_familial_controls_with_variants;
	int _show_external_controls_with_variants;
	bool _show_all_variants;
	bool _use_genotyper_reads;
//	bool _lenient_parsing;

	/* input-output file options */
	vector<string> _vcf_input_files;
	string _pedigree_input_file;
	string _output_file;
	vector<string> *_extra_info_fields;
	vector<string> *_extra_format_fields;
	vector<string> *_required_info_fields;
	vector<string> *_required_format_fields;
//	vector<string> *_format_filters;
	map<string, string> *_allInfoDescriptions;
	map<string, string> *_allFormatDescriptions;

	/* Boost options variable */
	options_description * _all;

	/* Singleton stuff */
	static CRunSettings* _instance;
	bool _isInitialized;

	/**
	 * (EMPTY) CONSTRUCTOR
	 */
	CRunSettings() {
		_all = new options_description(VERSION);
		_isInitialized = false;
		_extra_info_fields = new vector<string>();
		_extra_format_fields = new vector<string>();
		_required_info_fields = new vector<string>();
		_required_format_fields = new vector<string>();
//		_format_filters = new vector<string>();
		_allInfoDescriptions = new map<string, string>();
		_allFormatDescriptions = new map<string, string>();
	}

	~CRunSettings();

};

#endif /* CRUN_SETTINGS_H_ */
