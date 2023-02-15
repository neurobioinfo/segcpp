/*
 * ITypes.h
 *
 *  Created on: 2011-11-09
 *  Author: Alexandre Dionne-Laporte
 */

#ifndef ITYPES_H_
#define ITYPES_H_

#include <string>

using namespace std;

enum VariantFieldType {
	FIELD_FLOAT = 0
	, FIELD_INTEGER = 1
	, FIELD_FLAG = 2
	, FIELD_STRING = 3
	, FIELD_UNKNOWN = -1
};

enum VariantFieldNumber {
	ALLELE_NUMBER = -2
	, GENOTYPE_NUMBER = -1
	, UNDEFINED = -3
};

enum VariantZygosity {
	HOMOZYGOUS = 1
	, HETEROZYGOUS = 2
	, TRIZYGOUS = 3
	, WILDTYPE = 100
	, WILDTYPE_FILTERED = -100
	, HOMOZYGOUS_FILTERED = -1
	, HETEROZYGOUS_FILTERED = -2
	, TRIZYGOUS_FILTERED = -3
	, ZYGOSITY_UNKNOWN = 0
};

enum VariantType {
	UNKNOWN = 0,	/* huh ?! */
	SNV = 1,		/* snip */
	INDEL = 2,		/* insertion / deletion */
	CNV = 3,		/* copy number variant */
	SUBST = 4,		/* substitution */
	SV = 5,			/* structural variant */
	TRANS = 6		/* translocation */
};

const char output_seperator[51] = "--------------------------------------------------";

#endif /* ITYPES_H_ */

