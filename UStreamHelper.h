/**
 * Taken mostly from 'vcflib'
 */

#ifndef USTREAM_HELPER_
#define USTREAM_HELPER_

#include "ITypes.h"

#include<iterator>

template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::map<std::string, string>& m) {
	for ( std::map<string, string>::const_iterator it = m.begin(); it != m.end(); it++ ) {
		os << it->first << " " << it->second << endl;
	}
	return os;
}

std::ostream& operator<<(std::ostream& os, const std::map<std::string, VariantFieldType>& m)
{
	for ( std::map<string, VariantFieldType>::const_iterator it = m.begin(); it != m.end(); it++ ) {
		os << it->first << " " ;
	    switch (it->second) {
	        case FIELD_INTEGER:
	            os << "integer";
	            break;
	        case FIELD_FLOAT:
	        	os << "float";
	            break;
	        case FIELD_BOOL:
	        	os << "bool";
	            break;
	        case FIELD_STRING:
	        	os << "string";
	            break;
	        default:
	        	os << "unknown";
	            break;
	    }
	    os << endl;
	}
    return os;
}

ostream& operator<<(ostream& out, VariantFieldType type) {
    switch (type) {
        case FIELD_INTEGER:
            out << "integer";
            break;
        case FIELD_FLOAT:
            out << "float";
            break;
        case FIELD_BOOL:
            out << "bool";
            break;
        case FIELD_STRING:
            out << "string";
            break;
        default:
            out << "unknown";
            break;
    }
    return out;
}

ostream& operator<<(ostream& out, VariantFieldNumber number) {
    switch (number) {
        case ALLELE_NUMBER:
            out << "\"allele number\"";
            break;
        case GENOTYPE_NUMBER:
            out << "\"genotype number\"";
            break;
        case UNDEFINED:
            out << "undefined";
            break;
        default:
            out << (int)number;
            break;
    }
    return out;
}


#endif /* USTREAM_HELPER_ */
