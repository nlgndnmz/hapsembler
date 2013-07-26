/*****************************************************************************

	Part of Hapsembler package. See the README file for more information.
    Copyright (C) 2011-2013,  Nilgun Donmez <nild@cs.toronto.edu>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 ****************************************************************************/

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <climits>
#include <iostream>

#include "HapSuite.h"
#include "Scaffer.h"

using namespace std;

void usage(const char * progName, const char * badArg, const char * msg)
{
	cerr << endl << badArg << " " << msg << endl;
	cerr << "Please see the README file for more information. Version " << SCARPAVERSION << "(hsv:" << HAPVERSION << ")" << endl << endl;

	cerr << "    " << progName << " -c <file> -l <file> -i <file> -o <file> \n\n\
OPTIONS \n\
    --contigs|-c contigs_filename \n\
        Fasta formatted file containing the contigs (required) \n\
\n\
    --library|-l libraries_filename \n\
        File containing the paired read library information (required) \n\
\n\
    --mappings|-i mappings_filename \n\
        File containing the read mappings (required) \n\
\n\
    --output|-o output_filename \n\
        Output file for scaffolds (required) \n\
\n\
    --calibrate|-b [yes|no] \n\
        If set to yes, re-calculates the mean and standard deviation of each library. Default is yes. \n\
\n\
    --min_support N \n\
        Sets the minimum number of mate links to connect two contigs to N. Default value is 2. \n\
\n\
    --max_removal N \n\
        Sets the maximum number of contigs that can be removed during the orientation step to N. Default value is 6. \n\
\n\
    --min_contig N \n\
        Sets the minimum contig size (in bp) to be used in scaffolding to N. Default value is 100bp. \n\n";

    exit(0);
}

int main(int argc, char **argv)
{
	if(argc < 3)
		usage(argv[0], " ", " ");

	int minSupport = 2;
	int maxOverlap = 50;
	int maxDegree = 50;
	int minContig = 100;
	int maxK = 6;

	bool calibrate = true;

	char * mappingsFilename = NULL;
	char * contigsFilename = NULL;
	char * outputFilename = NULL;
	char * libraryFilename = NULL;

	for(int i=1; i<(argc-1); i+=2)
	{
		if(strcmp(argv[i], "--mappings") == 0 || strcmp(argv[i], "-i") == 0)
			mappingsFilename = argv[i+1];

		else if(strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "-o") == 0)
			outputFilename = argv[i+1];

		else if(strcmp(argv[i], "--contigs") == 0 || strcmp(argv[i], "-c") == 0)
			contigsFilename = argv[i+1];

		else if(strcmp(argv[i], "--library") == 0 || strcmp(argv[i], "-l") == 0)
			libraryFilename = argv[i+1];

		else if(strcmp(argv[i], "--calibrate") == 0 || strcmp(argv[i], "-b") == 0)
		{
			if(strcmp(argv[i+1], "no") == 0)
				calibrate = false;
		}

		else if(strcmp(argv[i], "--max_removal") == 0)
			sscanf( argv[i+1], "%d", &maxK);

		else if(strcmp(argv[i], "--min_contig") == 0)
			sscanf( argv[i+1], "%d", &minContig);

		else if(strcmp(argv[i], "--min_support") == 0)
			sscanf( argv[i+1], "%d", &minSupport);

		else
			usage(argv[0], argv[i], "unknown option. Please check your command.");
	}

	if(INT_MAX < 2097152000 || (LONG_MAX-10) < INT_MAX )	
		usage(argv[0], " ", "Integral types are too small! Please re-compile using a more recent compiler.");

	time_t rawtime;
	struct tm * timeinfo;

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	clock_t start, stop;
	start = clock();

	cout << endl << argv[0] << " started at " << asctime(timeinfo) << endl << "Arguments given: ";
	for(int i=1; i<argc; i++)
		cout << argv[i] << " ";
	cout << endl << endl;

	try
	{
		Scaffer * scaff = new Scaffer(maxDegree, maxOverlap, minContig, maxK);
		if(libraryFilename == NULL) usage(argv[0], " ", "Error: No library file is given!");
		if(contigsFilename == NULL) usage(argv[0], " ", "Error: No contigs file is given!");
		if(mappingsFilename == NULL) usage(argv[0], " ", "Error: No mappings file is given!");
		if(outputFilename == NULL) usage(argv[0], " ", "Error: No output file is given!");

		scaff->read_info(libraryFilename, minSupport);
		long int totalLength = scaff->read_mappings(mappingsFilename, outputFilename, calibrate);
		scaff->write_scaffolds(contigsFilename, outputFilename, totalLength);

		delete scaff;
	}
	catch(const char * str)
	{
		cerr << "Exception caught in " << argv[0] <<" : "<< str << endl;
	}

	time_t rawtime2;
	struct tm * timeinfo2;

	time ( &rawtime2 );
	timeinfo2 = localtime ( &rawtime2 );

	/* Stop timer */
	stop = clock();
	double t = (double) (stop-start)/CLOCKS_PER_SEC;

	cout << endl << argv[0] << " completed at " << asctime(timeinfo2) << " ( " << t << " )" << endl << endl;

	return 0;
}
