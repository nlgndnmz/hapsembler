/*****************************************************************************
    $Author: Nilgun Donmez $
    $Date: 2012-10-07 18:54:28 -0400 (Sun, 07 Oct 2012) $

	Part of Hapsembler package. See the README file for more information.
    Copyright (C) 2011,  Nilgun Donmez <nild@cs.toronto.edu>

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
#include <ctime>
#include <iostream>
#include <fstream>
#include <climits>

#include "HapSuite.h"
#include "HapUtils.h"
#include "OvlGraph.h"

using namespace std;

//***************************************************************
//    Methods
//***************************************************************

void usage(char * progname)
{
	cerr << "Please see the README file for more information. Version " << HAPVERSION << endl << endl;
    cerr << "\n\
SYNOPSIS \n\
    hapsemblr -r <file> -c <file> -g <genome size>\n\
\n\
OPTIONS \n\
    --prefix|-r reads_filename \n\
        Prefix of the files produced by overlappr (required) \n\
\n\
    --contigs|-c contigs_filename \n\
        Output filename for contigs (required) \n\
\n\
    --genome|-g genome_size \n\
        Estimated genome size in KILO base pairs (required) \n\
\n\
    --library|-l library_filename \n\
        File containing information about libraries. If unset, all reads are taken to be single production. \n\
\n\
    --onestrand|-a [yes|no] \n\
        If set to yes, performs strand specific assembly. This option can not be turned on if -l option is set. Default value is no. \n\
\n\
    --calibrate|-b [yes|no] \n\
        If set to yes, re-calculates the mean and standard deviation of each library. The default is yes. \n\n";

    exit(0);
}

int get_reads(OvlGraph *& OG, char * prefix, char * infoFilename, bool onestrand)
{
	string s(prefix);
	s += ".reads";

	char * readFilename = new char[s.size()+1];
	strcpy(readFilename, s.c_str());

	ifstream readFile;
	open_n_check(readFile, readFilename);

	int read_id = 1;
	int read_len = 0;
	int max_read_size = 0;
	int num_reads = 0;

	if( !(readFile >> num_reads) )
		throw "Can not read the reads file!";

	OG = new OvlGraph(num_reads, onestrand);

	while( readFile >> read_len )
	{
		OG->add_node(read_id++, read_len);
		if(read_len > max_read_size)
			max_read_size = read_len;
	}
	check_n_close(readFile);
	delete [] readFilename;

	if(infoFilename != NULL)
	{
		int num_pairs = 0;
		int num_libs = 0;

		ifstream infoFile;

		int start, end, insert_size, deviation, ori;

		open_n_check(infoFile, infoFilename);
		while( infoFile >> start >> end >> insert_size >> deviation >> ori )
		{
			if(insert_size < 16000)	// we can not process anything larger
				num_libs++;
			if( (end - start + 1)%2 == 1 )
				throw "Library contains odd number of reads!";
			num_pairs += (end - start + 1)/2;
		}
		check_n_close(infoFile);

		OG->set_pairs(num_pairs, num_libs);

		num_libs = 0;
		ifstream infoFile2;
		open_n_check(infoFile2, infoFilename);
		while( infoFile2 >> start >> end >> insert_size >> deviation >> ori )
		{
			int L = 0;
			if(insert_size < 16000)
			{
				OG->Libs[++num_libs].insertSize = insert_size;
				OG->Libs[num_libs].deviation = (short int) deviation;
				OG->Libs[num_libs].sdwidth = 3;
				OG->Libs[num_libs].maxDist = insert_size + 3*deviation;  // temporary value
				OG->Libs[num_libs].coverage = 1.0;
				OG->Libs[num_libs].orient = 1;		// field reserved for future

				L = num_libs;
			}

			for(int i=start+1; i<end; i+=2)
				OG->add_pair(i/2, L);
		}
		check_n_close(infoFile2);
	}

	return max_read_size;
}

void get_overlaps(OvlGraph * OG, char * prefix, int & min_overlap)
{
	string s(prefix);
	s += ".overlaps";

	char * ovlFilename = new char[s.size()+1];
	strcpy(ovlFilename, s.c_str());

	ifstream ovlFile;
	open_n_check(ovlFile, ovlFilename);

	int read1, read2;
    int comp2, st1, st2, e1, e2, diff, indels;

    double flexi = 1.0001;

    while( ovlFile >> read1 >> read2 >> comp2 >> st1 >> st2 >> e1 >> e2 >> diff >> indels )
    {
    	int overlap_len = int( (e1 - st1 + e2 - st2) / 2.0 );

    	if(((double)indels/overlap_len) + 1.0 > flexi)
			flexi = ((double)indels/overlap_len) + 1.0;

    	if(read1 < read2)
    	{
    		OG->add_arc(read1, read2, st1, st2, e1, e2, comp2, true);

    		if(overlap_len < min_overlap)
				min_overlap = overlap_len;
    	}
    }
	check_n_close(ovlFile);
	delete [] ovlFilename;

    OG->set_flexi(flexi);
}

//***************************************************************
//    Main
//***************************************************************
int main(int argc, char **argv)
{
	if(argc < 3)
	{    usage(argv[0]);    }

	bool onestrand = false;
	bool calibrate = true;
	int genome = 0;

	char * prefix = NULL;
	char * infoFilename = NULL;
	char * contigFilename = NULL;

	for(int i=1; i<(argc-1); i+=2)
	{
		if(strcmp(argv[i], "--prefix") == 0 || strcmp(argv[i], "-r") == 0)
			prefix = argv[i+1];

		else if(strcmp(argv[i], "--library") == 0 || strcmp(argv[i], "-l") == 0)
			infoFilename = argv[i+1];

		else if(strcmp(argv[i], "--contigs") == 0 || strcmp(argv[i], "-c") == 0)
			contigFilename = argv[i+1];

		else if(strcmp(argv[i], "--genome") == 0 || strcmp(argv[i], "-g") == 0)
			sscanf(argv[i+1], "%d", &genome);

		else if(strcmp(argv[i], "--calibrate") == 0 || strcmp(argv[i], "-b") == 0)
		{
			if(strcmp(argv[i+1], "no") == 0)
				calibrate = false;
		}
		else if(strcmp(argv[i], "--onestrand") == 0 || strcmp(argv[i], "-a") == 0)
		{
			if(strcmp(argv[i+1], "yes") == 0)
				onestrand = true;
		}
		else
		{
			cerr << "Unknown parameter! Please check your command." << endl;
			usage(argv[0]);
		}
	}

	if(UINT_MAX <= 2147483647)
	{
		cerr << "Integer types are too small! Please re-compile Hapsembler using a more recent compiler!" << endl;
		exit(0);
	}

	time_t rawtime;
	struct tm * timeinfo;

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	clock_t start, stop;
	start = clock();

	cout << argv[0] << " started at " << asctime(timeinfo) << endl << "Arguments: ";
	for(int i=1; i<argc; i++)
		cout << argv[i] <<" ";
	cout << endl << endl;

	if(infoFilename != NULL && onestrand)
	{
		cerr << "Strand specific assembly can not be turned on in paired mode!" << endl;
		usage(argv[0]);
	}

	if(prefix == NULL) { cerr << "Prefix for input files is not given!" << endl; usage(argv[0]); }
	if(contigFilename == NULL) { cerr << "Output filename is not given!" << endl; usage(argv[0]); }
	if(genome < 1) { cerr << "Genome size is too small or not given!" << endl; usage(argv[0]); }

	genome *= 1000;		// convert it to base pairs

	OvlGraph * OG = NULL;

	try
	{
		int min_overlap = get_reads(OG, prefix, infoFilename, onestrand);

		cout << "Finished reading the reads" << endl;

		get_overlaps(OG, prefix, min_overlap);

		cout << endl << "Processing read graph" << endl << endl;

		OG->sort_read_edges();
		OG->check_read_graph();

		OG->remove_contained_reads_alt();
		OG->check_read_graph();

		int n = 0;
		while( OG->venom() > 0 ) n++;
		cout << "Finished " << n << " iterations of edge splitting" << endl;
		OG->check_read_graph();

		OG->reduce_read_graph();
		OG->check_read_graph();

		while( OG->remove_read_buds() > 0 );
		OG->check_read_graph();

		cout << "Number of read bubbles removed: " << OG->remove_read_bubbles() << endl;
		int numActive = OG->check_read_graph();
		OG->read_chain_collapse();

		if(infoFilename == NULL)
		{
			OG->simplify_read_chains(max(2, int(100.0/(genome/(double)numActive))));
			OG->print_read_contigs(contigFilename);
		}
		else
		{
			if(calibrate)
				OG->calibrate_libraries(contigFilename);

			OG->build_pair_graph();

			cout << endl << "Processing mate graph" << endl << endl;

			OG->sort_mate_edges();
			OG->check_mate_graph();

			OG->remove_contained_mates();
			OG->check_mate_graph();

			OG->reduce_mate_graph();
			OG->check_mate_graph();

			while( OG->remove_mate_buds() > 0 );
			OG->check_mate_graph();

			cout << "Number of mate bubbles removed: " << OG->remove_mate_bubbles() << endl;
			numActive = OG->check_mate_graph();

			OG->bridge_matepairs();
			OG->check_mate_graph();

			OG->mate_chain_collapse();
			OG->simplify_mate_chains(max(2, int(200.0/(genome/(double)numActive))));

			cout << "Number of gaps: " << OG->print_mate_contigs(contigFilename) << endl;
		}
	}
	catch(const char * str)
	{
		fprintf(stderr, "Exception caught: %s", str);
	}

	time_t rawtime2;
	struct tm * timeinfo2;

	time ( &rawtime2 );
	timeinfo2 = localtime ( &rawtime2 );

	stop = clock();
	double t = (double) (stop-start)/CLOCKS_PER_SEC;

	delete OG;

	cout << endl << argv[0] << " completed at " << asctime(timeinfo2) << " ( " << t << " )" << endl << endl;

    return 0;
}
