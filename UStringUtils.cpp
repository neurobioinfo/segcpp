/**
 * Taken mostly from 'vcflib'
 */

#include "UStringUtils.h"

#include <string.h>
#include <iostream>

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

std::vector<std::string> &split(const std::string &s, const std::string &delims, std::vector<std::string> &elems) {
    char* tok;
    char cchars [s.size()+1];
    char* cstr = &cchars[0];
    strcpy(cstr, s.c_str());
    tok = strtok(cstr, delims.c_str());
    while (tok != NULL) {
        elems.push_back(tok);
        tok = strtok(NULL, delims.c_str());
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, const std::string &delims) {
    std::vector<std::string> elems;
    return split(s, delims, elems);
}

std::vector<std::string> splitWithEmpty(const std::string &record, const std::string &token) {
    vector<string> results;
    size_t startPos = 0;
    size_t pos = 0;

    // Step: If either argument is empty then return an empty vector.
    if (record.length() == 0 || token.length() == 0) {
        return results;
    }

    // Step: Go through the record and split up the data.
    while(startPos < record.length()) {
        pos = record.find(token, startPos);
        if (pos == string::npos) {
            break;
        }

        results.push_back(record.substr(startPos, pos - startPos));
        startPos = pos + token.length();
    }

    // Step: Get the last (or only bit).
    results.push_back(record.substr(startPos, record.length() - startPos));

    // Step: Return the results of the split.
    return results;
}

std::vector<std::string> splitStrDelim(const std::string &s, const std::string& delimiter) {
	std::vector<std::string> elems;
	size_t pos = 0;
	std::string token;
	std::string s_copy = s;
	if( s_copy.find(delimiter) == std::string::npos ) { elems.push_back(s); }
	else {
		while ((pos = s_copy.find(delimiter)) != std::string::npos) {
			token = s_copy.substr(0, pos);
			elems.push_back(token);
			s_copy.erase(0, pos + delimiter.length());
		}
	}
    return elems;
}

VariantFieldType str2VariantFieldType(const std::string &typeStr) {
	if (typeStr == "Integer") {
		return FIELD_INTEGER;
	} else if (typeStr == "Float") {
		return FIELD_FLOAT;
	} else if (typeStr == "String") {
		return FIELD_STRING;
	} else if (typeStr == "Flag") {
		return FIELD_FLAG;
	} else {
		return FIELD_UNKNOWN;
	}
}
const std::string variantFieldType2str(const VariantFieldType &type) {
	if ( type == FIELD_INTEGER ) {
		return "Integer";
	}
	else if ( type == FIELD_FLOAT ) {
		return "Float";
	}
	else if ( type == FIELD_FLAG ) {
		return "Flag";
	}
	else if ( type == FIELD_STRING ) {
		return "String";
	}
	else {
		return "Uknown";
	}
}

int str2VariantFieldNumber(const std::string &numStr) {
	if (numStr == "A") {
		return ALLELE_NUMBER;
	} else if (numStr == "G") {
		return GENOTYPE_NUMBER;
	} else if (numStr == ".") {
		return UNDEFINED;
	} else {
		int retVal;
		convert(numStr, retVal);
		return retVal;
	}
}
const std::string variantFieldNumber2str(const int &num) {
	if ( num == ALLELE_NUMBER ) {
		return "A";
	}
	else if ( num == GENOTYPE_NUMBER ) {
		return "G";
	}
	else if ( num == UNDEFINED ) {
		return ".";
	}
	else {
		return convert(num);
	}
}

const std::string variantZygosity2str(const VariantZygosity &zygosity) {
	if ( zygosity == HOMOZYGOUS ) {
		return "Homozygote";
	}   else if ( zygosity == HETEROZYGOUS ) {
		return "Heterozygote";
	}	else if ( zygosity == TRIZYGOUS ) {
		return "Trizygote";
	}	else if ( zygosity == WILDTYPE ) {
		return "Wildtype";
	}	else if ( zygosity == WILDTYPE_FILTERED ) {
		return "Filtered Wildtype";
	}	else if ( zygosity == HOMOZYGOUS_FILTERED ) {
		return "Filtered Homozygote";
	}	else if ( zygosity == HETEROZYGOUS_FILTERED ) {
		return "Filtered Heterozygote";
	}	else if ( zygosity == TRIZYGOUS_FILTERED ) {
		return "Filtered Trizygote";
	}	else {
		return "Uknown";
	}
}
