/**
 * Taken mostly from 'vcflib'
 */

#ifndef USTRING_UTILS_
#define USTRING_UTILS_

#include "ITypes.h"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

/**
 * converts the string into the specified type,
 * setting r to the converted value and returning
 * true/false on success or failure
 */
template<typename T>
bool convert(const std::string& s, T& r) {
    std::istringstream iss(s);
    iss >> r;
    return (iss.fail() || ((std::size_t) iss.tellg()) != s.size()) ? false : true;
}

template<typename T>
std::string convert(const T& r) {
    std::ostringstream iss;
    iss << r;
    return iss.str();
}

VariantFieldType str2VariantFieldType(const std::string &typeStr);
const std::string variantFieldType2str(const VariantFieldType &type);

int str2VariantFieldNumber(const std::string &numStr);
const std::string variantFieldNumber2str(const int &num);

//VariantZygosity str2VariantZygosity(const std::string &zygosityStr);
const std::string variantZygosity2str(const VariantZygosity &zygosity);

/**
 * functions to split a string by a specific delimiter
 */

/**
 * join a vector of elements by a delimiter object.  ostream<< must be defined
 * for both class S and T and an ostream, as it is e.g. in the case of strings
 * and character arrays
 */

template<class S, class T>
std::string join(std::vector<T>& elems, S& delim) {
    std::stringstream ss;
    typename std::vector<T>::iterator e = elems.begin();
    ss << *e++;
    for (; e != elems.end(); ++e) {
        ss << delim << *e;
    }
    return ss.str();
}

// thanks to Evan Teran, http://stackoverflow.com/questions/236129/how-to-split-a-string/236803#236803

/**
 * split a string on a single delimiter character (delim)
 */
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);

/**
 * split a string on any character found in the string of delimiters (delims)
 */
std::vector<std::string> &split(const std::string &s, const std::string& delims, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, const std::string& delims);
std::vector<std::string> splitWithEmpty(const std::string &record, const std::string &token);

/**
 * split a string on all characters found in the string of delimiters (delims)
 */
std::vector<std::string> splitStrDelim(const std::string &s, const std::string& delims);

#endif /* USTRING_UTILS_ */
