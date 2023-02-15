/*
 * CVariantInSample.cpp
 *
 *  Created on: 2011-02-15
 *  Author: Alexandre Dionne-Laporte
 */

#include "CVariantInSample.h"
#include "CSample.h"
#include "CVariant.h"
#include "CRunSettings.h"

#include "iostream"

CVariantInSample::CVariantInSample(
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
) :
_sample(sample),
_variant(variant),
_genotype_quality(genotype_quality),
_genotype(zygosity),
_depth(depth),
_genotyper_reads_ref(genotyper_reads_ref),
_genotyper_reads_var(genotyper_reads_var),
_total_reads_ref_plus(total_reads_ref_plus),
_total_reads_ref_minus(total_reads_ref_minus),
_total_reads_var_plus(total_reads_var_plus),
_total_reads_var_minus(total_reads_var_minus)

{
	_sample->addVariantInSample(this);
	_variant->addVariantInSample(this);
	_formats = new map<string, string >();
}

CVariantInSample::~CVariantInSample() {
	delete _sample;
	delete _variant;
	_formats->clear();
	delete _formats;
//	if (CRunSettings::getInstance()->get_verbose_level() > CRunSettings::VERBOSE_LOG)
//		cerr << "CVariantInSample::~CVariantInSample" << endl;
}
