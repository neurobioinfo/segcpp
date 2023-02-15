/*
 * CSample.h
 *
 * Created on: 2011-02-01
 * Author: Alexandre Dionne-Laporte
 */

#ifndef CSAMPLE_H_
#define CSAMPLE_H_

class CFamily;
class CVariantInSample;

#include <map>
#include <vector>
#include <string>
#include <iosfwd>

using namespace std;

class CSample {
	friend std::ostream & operator << (std::ostream &, const CSample &);

protected:
	string _name;
	bool _isProband;
	CFamily * _family;
	vector<CVariantInSample*> * _variantsInSample;

private:
	static map<string, CSample*> * _samples;

public:
	//STATIC
	static CSample* getSample(string);
	static vector<CSample*> getSamples();

	//CONSTRUCTORS
	CSample(string name, bool isProband, CFamily * family);
	~CSample();

	//GETTERS
	inline string getName () const { return _name; };
	inline bool isProband () const { return _isProband; };
	inline CFamily* getFamily () const { return _family; }
	inline vector<CVariantInSample*> * getVariantsInSample() const { return _variantsInSample; }

	//SETTERS
	inline void addVariantInSample (CVariantInSample* variantInSample) { _variantsInSample->insert(_variantsInSample->end(), variantInSample); };

};

#endif /* CSAMPLE_H_ */
