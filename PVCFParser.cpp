#include <iostream>
#include <iterator>
#include <vector>
#include <set>
#include <fstream>

#include "CVCFFile.h"
#include "CRunSettings.h"
#include "UStringUtils.h"
#include "UPedFileParser.h"

//Boost stuff
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace boost;
using namespace boost::program_options;

using namespace std;

int main( int ac, char *av[ ] )
{
//	string tata("aaatata");
//	vector<string> v = split(tata, 't');
//	cout << v << endl;


	ifstream toto(av[2]);
	vector<CSample*> * samples = UPedFileParser::parseFile(toto);
	CVCFFile vcfFile (samples);

	vcfFile.openFile(string(av[1]));
//	cout << "FormatTypes" << endl << vcfFile.getFormatTypes() << endl;
//	cout << "FormatDescriptions" << endl << vcfFile.getFormatDescriptions() << endl;
//	cout << "InfoTypes" << endl << vcfFile.getInfoTypes() << endl;
//	cout << "InfoDescriptions" << endl << vcfFile.getInfoDescriptions() << endl;
//	cout << "Samples: " << endl << vcfFile.getSampleNames() << endl;

	vector<CVariant*>* variants = vcfFile.parseFile();

	for ( vector<CVariant*>::iterator i = variants->begin() ; i != variants->end(); i++) {
		cout << **i << endl;
		cout << "=================" << endl;
	}

	return 1;
}


