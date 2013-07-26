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

#include "HapSuite.h"
#include "Dna.h"

using namespace std;

void usage(const char * progName, const char * badArg, const char * msg)
{
	cerr << endl << badArg << " " << msg << endl;
#ifdef USEOPENMP
	cerr << "Please see the README file for more information. Version " << ENCOREVERSION << "(hsv:" << HAPVERSION << ") (OpenMP support: yes) " << endl << endl;
#else
	cerr << "Please see the README file for more information. Version " << ENCOREVERSION << "(hsv:" << HAPVERSION << ") (OpenMP support: no) " << endl << endl;
#endif

#ifdef PREPROCR
	cerr << "    " << progName << " -p <platform> -f <file> -o <file> " << endl << "\n\
OPTIONS \n\
    --platform|-p [illumina|fourfivefour] \n\
        Define the type of the platform the reads are produced from (required) \n\
\n\
    --fastq|-f fastq_filename \n\
        Fastq formatted input file (required) \n\
\n\
    --fastq2|-x fastq_filename \n\
        The second reads in fastq format for paired libraries. The order of the reads should match the order of the reads in the file given with --fastq option. \n\
\n\
    --output|-o output_filename \n\
        Output filename. Default is standard output. \n\
\n\
    --revcomp|-n [1|2|3] \n\
        Reverse complement the first (1), the second (2) or both (3) reads in a read pair. If this option is set all reads are taken to be paired. \n\
\n\
    --phred|-d N \n\
        Set the phred offset for the quality values to N. Default value is 33 for fourfivefour and 64 for illumina. \n\
\n\
    --threshold|-s E \n\
        Set the threshold for trimming. E must be a real number between 0.0 (no trimming) and 0.5 (most trimming). Default values are 0.05 for illumina and 0.1 for fourfivefour. \n\n";
#endif

#ifdef CORRECTR
	cerr << "    " << progName << " -p <platform> -f <file> -o <file> -g <genome> " << endl << "\n\
OPTIONS \n\
    --platform|-p [illumina|fourfivefour] \n\
        Define the type of the platform the reads are produced from (required) \n\
\n\
    --fastq|-f fastq_filename \n\
        Fastq formatted input file (required) \n\
\n\
    --output|-o output_filename \n\
        Output filename for corrected reads (required) \n\
\n\
    --genome|-g V \n\
        Estimated genome size in KILO base pairs. (required) \n\
\n\
    --nthreads|-t K \n\
        Use K number of threads (ignored if program is compiled without OpenMP) \n\
\n\
    --onestrand|-a [yes|no] \n\
        If set to yes, the reads are treated as single stranded. Default value is no. \n\
\n\
    --epsilon|-e X (real number) \n\
        Set the expected discrepancy (mismatches+indels) rate. X must be a real number between 0.01 and 0.09. Default values for illumina and fourfivefour are 0.04 and 0.06 respectively. \n\
\n\
    --phred|-d N \n\
        Set the phred offset for the quality values to N. Default value is 33 for fourfivefour and 64 for illumina. \n\n";
#endif

#ifdef OVERLAPPR
	cerr << "    " << progName <<  " -p <platform> -f <file> -o <prefix> -g <genome>" << endl << "\n\
OPTIONS \n\
    --platform|-p [illumina|fourfivefour] \n\
        Define the type of the platform the reads are produced from (required) \n\
\n\
    --fastq|-f fastq_filename \n\
        Fastq formatted input file (required) \n\
\n\
    --output|-o prefix \n\
        Prefix for output files \n\
\n\
    --genome|-g V \n\
        Estimated genome size in KILO base pairs. (required) \n\
\n\
    --nthreads|-t K \n\
        Use K number of threads (ignored if program is compiled without OpenMP) \n\
\n\
    --onestrand|-a [yes|no] \n\
        If set to yes, the reads are treated as single stranded. Default value is no. \n\
\n\
    --epsilon|-e X (real number) \n\
        Set the expected discrepancy (mismatches+indels) rate. X must be a real number between 0.01 and 0.09. Default values for illumina and fourfivefour are 0.04 and 0.06 respectively. \n\n";
#endif

#ifdef CONSENSR
	cerr << "    " << progName << " -p <platform> -f <file> -c <file> -o <file> " << endl << "\n\
OPTIONS \n\
    --platform|-p [illumina|fourfivefour] \n\
        Define the type of the platform the reads are produced from (required) \n\
\n\
    --fastq|-f fastq_filename \n\
        Fastq formatted input file (required) \n\
\n\
    --contigs|-c contigs_filename \n\
        Contig file produced by hapsemblr (required) \n\
\n\
    --output|-o output_filename \n\
        Output file for contigs (required) \n\
\n\
    --min-size|-m N (integer) \n\
        Set the minimum size of a contig to be reported to N (in bp). Default value is 200. \n\
\n\
    --phred|-d N \n\
        Set the phred offset for the quality values to N. Default value is 33 for fourfivefour and 64 for illumina. \n\n";
#endif

    exit(0);
}

int main(int argc, char **argv)
{
	if(argc < 3)
		usage(argv[0], " ", " ");

	// misc parameters
	int kmer_size = 14;
	int min_contig_size = 200;
	double trim_threshold = -0.1;

	// default Smith-Waterman parameters
	int match_score = 1;
	int mismatch_score = -2;
	int indel_score = -6;

	// platform specific defaults (will be overridden)
	double epsilon = -1.0;
	int phred_offset = 0;
	int gap_quality = -1;

	int max_read_size = 32000;

	// default Naive Bayes arguments
	int max_quality = 100;
	int num_nucleotides = 4;
	int max_steps = 1000;
	double prior = 0.999;

	long int genome = 0;
	Sequencer platform = NOPLATFORM;

	char * fastqFilename = NULL;
	char * pairedFilename = NULL;
	char * contigsFilename = NULL;
	char * outputFilename = NULL;

	int nthreads = 1;
	bool onestrand = false;
	int revcomp = 0;

	for(int i=1; i<(argc-1); i+=2)
	{
		if(strcmp(argv[i], "--fastq") == 0 || strcmp(argv[i], "-f") == 0)
			fastqFilename = argv[i+1];

		else if(strcmp(argv[i], "--fastq2") == 0 || strcmp(argv[i], "-x") == 0)
			pairedFilename = argv[i+1];

		else if(strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "-o") == 0)
			outputFilename = argv[i+1];

		else if(strcmp(argv[i], "--contigs") == 0 || strcmp(argv[i], "-c") == 0)
			contigsFilename = argv[i+1];

		else if(strcmp(argv[i], "--genome") == 0 || strcmp(argv[i], "-g") == 0)
			sscanf( argv[i+1], "%ld", &genome);

		else if(strcmp(argv[i], "--revcomp") == 0 || strcmp(argv[i], "-n") == 0)
			sscanf( argv[i+1], "%d", &revcomp);

		else if(strcmp(argv[i], "--nthreads") == 0 || strcmp(argv[i], "-t") == 0)
			sscanf( argv[i+1], "%d", &nthreads);

		else if(strcmp(argv[i], "--platform") == 0 || strcmp(argv[i], "-p") == 0)
		{
			if(strcmp(argv[i+1], "illumina") == 0)
				platform = ILLUMINA;
			else if(strcmp(argv[i+1], "fourfivefour") == 0)
				platform = FOURFIVEFOUR;
			else
				usage(argv[0], argv[i+1], "unknown option. Please check your command.");
		}
		else if(strcmp(argv[i], "--onestrand") == 0 || strcmp(argv[i], "-a") == 0)
		{
			if(strcmp(argv[i+1], "yes") == 0)
				onestrand = true;
		}

		else if(strcmp(argv[i], "--epsilon") == 0 || strcmp(argv[i], "-e") == 0)
			sscanf( argv[i+1], "%lf", &epsilon);

		else if(strcmp(argv[i], "--phred") == 0 || strcmp(argv[i], "-d") == 0)
			sscanf( argv[i+1], "%d", &phred_offset);

		else if(strcmp(argv[i], "--min-size") == 0 || strcmp(argv[i], "-m") == 0)
			sscanf( argv[i+1], "%d", &min_contig_size);

		else if(strcmp(argv[i], "--threshold") == 0 || strcmp(argv[i], "-s") == 0)
			sscanf( argv[i+1], "%lf", &trim_threshold);

		else
			usage(argv[0], argv[i], "unknown option. Please check your command.");
	}

	switch(platform)		// set up the platform specific parameters
	{
		case ILLUMINA:
			gap_quality = 70;
			if(epsilon < 0.0) epsilon = 0.04;
			if(phred_offset == 0) phred_offset = 64;
			if(trim_threshold < 0.0) trim_threshold = 0.05;
			break;

		case FOURFIVEFOUR:
			gap_quality = 10;
			if(epsilon < 0.0) epsilon = 0.06;
			if(phred_offset == 0) phred_offset = 33;
			if(trim_threshold < 0.0) trim_threshold = 0.1;
			break;

		case NOPLATFORM:
			usage(argv[0], " ", "Please specify a platform!");
	}

	// make sure the parameters are somewhat sane
	if(phred_offset < 33 || phred_offset > 80) usage(argv[0], " ", "Phred offset should be between 33 and 80");
	if(epsilon < 0.01) epsilon = 0.01;
	if(epsilon > 0.09) epsilon = 0.09;
	if(nthreads < 1) nthreads = 1;
	int thread_limit = DNA::alphabet_size * DNA::alphabet_size;
	if(nthreads > thread_limit) nthreads = thread_limit;
	if(trim_threshold > 0.5) trim_threshold = 0.5;
	if(revcomp < 0 || revcomp > 3) revcomp = 0;
	if(genome > 3200000) usage(argv[0], " ", "Genome size exceeds 3,200,000kbp. Please set the genome size in KILO base pairs.");

	if(INT_MAX < 2097152000 || (LONG_MAX-10) < INT_MAX )	
		usage(argv[0], " ", "Integral types are too small! Please re-compile using a more recent compiler.");

	time_t rawtime;
	struct tm * timeinfo;

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	clock_t start, stop;
	start = clock();

#ifndef USEOPENMP
	if(nthreads > 1)
		cerr << "You have requested multiple threads but the program is compiled without OpenMP support. Please recompile to enable multi-threading." << endl << endl;
#endif

	cerr << endl << argv[0] << " started at " << asctime(timeinfo) << endl << "Arguments given: ";
	for(int i=1; i<argc; i++)
		cerr << argv[i] << " ";
	cerr << endl << endl;

	if(fastqFilename == NULL) usage(argv[0], " ", "Error: No fastq file is given!");

#ifndef PREPROCR
	if(outputFilename == NULL) usage(argv[0], " ", "Error: No output file is given!");
#endif

	try
	{
		DNA * seqs = new DNA(kmer_size, max_read_size, phred_offset, platform);
		seqs->set_SmithWaterman(max_read_size, 1+2*kmer_size, match_score, mismatch_score,
			indel_score, gap_quality, epsilon, max_quality, num_nucleotides, max_steps, prior);

#ifdef PREPROCR
		seqs->trim_reads(fastqFilename, pairedFilename, outputFilename, kmer_size, trim_threshold, revcomp);
#endif

		genome *= 1000;		// convert it to base pairs

#ifdef CORRECTR
		seqs->isCorrectr = true;
		if(genome < 1000) usage(argv[0], " ", "Error: Genome size can not be less than 1kbp!");
		int coverage = seqs->read_reads(fastqFilename, genome, true, true);

		nthreads = seqs->overlap_reads(nthreads, outputFilename, onestrand, coverage);
		seqs->cat_reads(outputFilename, nthreads);
#endif

#ifdef OVERLAPPR
		seqs->isOverlappr = true;
		if(genome < 1000) usage(argv[0], " ", "Error: Genome size can not be less than 1kbp!");
		int coverage = seqs->read_reads(fastqFilename, genome, false, false, outputFilename);

		nthreads = seqs->overlap_reads(nthreads, outputFilename, onestrand, coverage);
		seqs->prune_overlaps(outputFilename, nthreads);
#endif

#ifdef CONSENSR
		seqs->isConsensr = true;
		if(contigsFilename == NULL) usage(argv[0], " ", "Error: No contigs file is given!");

		seqs->read_reads(fastqFilename, genome, true);
		seqs->write_contigs(contigsFilename, min_contig_size, outputFilename);
#endif

		delete seqs;
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

	cerr << endl << argv[0] << " completed at " << asctime(timeinfo2) << " ( " << t << " )" << endl << endl;

	return 0;
}
