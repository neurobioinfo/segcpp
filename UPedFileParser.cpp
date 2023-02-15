/*
 * CPedFile.cpp
 *
 *  Created on: 2011-10-31
 *  Author: Alexandre Dionne-Laporte
 */

#include "CFamily.h"
#include "UPedFileParser.h"
#include "CRunSettings.h"
#include "CSample.h"

#include <boost/algorithm/string.hpp>

using namespace std;

vector<CSample*> * UPedFileParser::parseFile (std::ifstream&  inStream) throw (EPedFileParser) {

/*
 * FORMAT:
 *

Family ID | Individual ID | Paternal ID | Maternal ID | Sex (1=male; 2=female; other=unknown) | Phenotype

FAM  ID			P		M		S P
Fam1 Sample1 	0		0		1 1
Fam1 Sample2 	0		0		2 1
Fam1 Sample3	Sample1	Sample2 1 1
Fam1 Sample4	Sample2	Sample2 1 2
Fam2 Sample5	0		0		1 1
Fam2 Sample6	Sample5	0 		2 1
Fam3 Sample7 	0 		0 		1 2
Fam3 Sample8 	0		0		1 2
Fam3 Sample9	0		0		1 2

 *
 */
	vector<CSample*> * samples = new vector<CSample*> ();
	string pedLine;

	while (inStream.good()) {
		getline(inStream, pedLine, '\n' );

		if( pedLine.empty() || pedLine.at(0) == '#') {
			continue;
		}

		std::vector<std::string> pedColumns;
		boost::split(pedColumns, pedLine, boost::is_any_of("\t "));

		if ( pedColumns.size() != 6 ) {
			throw EPedFileParser (string("Error in ped file, columns != 6" + pedLine));
		}

		CFamily* family = CFamily::getSetFamily( (string)(pedColumns.at(0)) );
		bool isProband = (atoi(((string)(pedColumns.at(5))).c_str()) == 2) ? true : false;
		CSample* sample = new CSample( (string)(pedColumns.at(1)), isProband, family );

		/* check that a sample is only included once in the ped file */
		for ( vector<CSample*>::iterator sampleIt = samples->begin(); sampleIt != samples->end(); sampleIt++ ) {
			CSample* currentSample = *sampleIt;
			if ( currentSample->getName().compare(sample->getName())==0) {
				//problem!
				throw EPedFileParser (string("Duplicate sample in ped file:" + sample->getName()));
			}

		}

		family->addSample(sample);
		if ( isProband ) { samples->insert ( samples->begin() , sample ); }
		else { samples->insert ( samples->end() , sample ); }

	}
	return samples;
}
