/*****************************************************************************
    $Author: Nilgun Donmez $
    $Date: 2012-09-02 18:28:18 -0400 (Sun, 02 Sep 2012) $

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

#ifndef HAP_UTILS_H
#define HAP_UTILS_H

#include <fstream>
#include <vector>

#include "GraphDef.h"

// opens an input file
void open_n_check(std::ifstream & fh, char * filename);

// opens an output file
void open_n_check(std::ofstream & fh, char * filename);

// closes an input file (ifstream)
void check_n_close(std::ifstream & fh);

// closes an output file (ofstream)
void check_n_close(std::ofstream & fh);

// swaps elements at indices qs1 and qs2 in que (also update parent pointers ptrs)
void swap_que_elts(int qs1, int qs2, QStruct * que, PStruct * ptrs);

// Required by qsort at various places, sorts in non-decreasing order
int compare(const void * a, const void * b);

// Used to report N50 and other assembly statistics by consensr and scaffr
int length_stats(int * lengths, int num, long int total, double frac=0.5);

// Required by the lexdfs method in OvlGraph.h
bool larger(const std::vector<int> & v1, const std::vector<int> & v2);

// Required by the lexdfs method in OvlGraph.h
bool larger_eq(const std::vector<int> & v1, const std::vector<int> & v2);

// used by correct_reads in Dna.h
void calculate_het_cutoff(int * table, int max_overlaps, double err, int max_t, double prob);

// used by OvlGraph.cpp and by Scaffer.cpp to adjust the library statistics
void bootstrap_library(long int * libMean, int * libSupport, int ** histogram,
 int histSize, Library * Libs, int k, bool upperOnly=false);

#endif // HAP_UTILS_H
