/*****************************************************************************
    $Author: Nilgun Donmez $
    $Date: 2012-08-29 21:24:20 -0400 (Wed, 29 Aug 2012) $

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

#ifndef DNA_H
#define DNA_H

#include <iostream>
#include <string>

#include "GraphDef.h"

class ContigNode;
class PathNode;
class DnaRead;
class KmerHash;
class Naive;
class ReadAln;
class ReadMatch;
class SmithWaterman;

class DNA
{
	public:

		DNA(int, int, int, Sequencer);
		void set_SmithWaterman(int, int, int, int, int, int, double, int, int, int, double);
		~DNA();

		void trim_reads(char *, char *, char *, int, double, int);

		int read_reads(char *, int, bool = false, bool = false, char * = NULL);
		int overlap_reads(int, char *, bool, int);
		void prune_overlaps(char *, int);
		void cat_reads(char *, int);

		void write_contigs(char *, int, char *);

		// macro mimickers
		bool isConsensr;
		bool isCorrectr;
		bool isOverlappr;

		const static int alphabet_size = 26;

	private:

		std::string read_next_fasta(std::ifstream &, DnaRead &, bool &);
		std::string read_next_fastq(std::ifstream &, DnaRead &, bool &, bool = true);

		bool handle_indels(int, ReadAln *, int, bool);

		bool * scan_contigs(char *, int);
		int write_contig(std::ofstream &, std::ofstream &, ContigTuple *, int, int);

		void correct_read(std::ofstream &, char *, int, ReadAln *, int, double *);
		int get_overlaps(int, int, ReadAln *, SmithWaterman *);

		void fill_in_kmers(int, int, DnaRead *, int);
		void count_qmers(DnaRead *, int, int);
		bool compare_reads(DnaRead &, bool, ReadMatch *, RMatch *, DnaRead *, SmithWaterman *);
		bool align_reads(DnaRead &, DnaRead &, bool, bool, double, ReadMatch *, SmithWaterman *);

		// the basics
		DnaRead * reads;
		char * deflines;
		char * sequences;
		char ** defs;

		int num_reads;
		int max_read_size;

		// alphabet related
		int * alphabet;
		char * alph;
		char * capit;

		// overlapper/mapper related fields
		KmerHash * KH;
		int kmer_size;

		// helper classes
		SmithWaterman * SW;
		Naive * NB;

		Sequencer platform;
		int expected_overlaps;
		int rmax_overlaps;
};

#endif // DNA_H
