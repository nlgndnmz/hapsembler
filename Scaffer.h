/*****************************************************************************
    $Author: Nilgun Donmez $
    $Date: 2011-12-06 17:19:54 -0500 (Tue, 06 Dec 2011) $

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

#ifndef SCAFFER_H
#define SCAFFER_H

#include <string>

#include "GraphDef.h"

class ContigNode;
class ContigEdge;

class Scaffer
{
	public:

		Scaffer(int, int, int, int);
		~Scaffer();

		void read_info(char *, int);
		int read_mappings(char *, char *, bool);
		int write_scaffolds(char *, char *, long int);

	private:

		void adjust_libraries(ReadTag *, bool, long int, int, ContigTag *, double);

		int fill_contigs(ReadTag *, int);
		void chain_collapse();
		void write_scafftmp(char *);

		int peekaboo(bool, int);
		int prune_graph();

		bool assign_components(int, int &, int=1);
		bool process_biconnected_components(int * que, int numMembers, int * map, int * map2, int step, int lib, int &);

		void find_articulation_points(int * que, int * map, int numMembers);
		void mark_articulation_edges(int ** edges, ContigEdge ** edgePtrs, int numMembers, int n);

		bool orient_component(int *, int, int, int, int *, int);
		void feedback_arc_set(int *, int, int);
		int solve_lp(int *, int, int *, int, int);

		int ** floyd_warshall(int * que, int n, int * map, int comp);

		ReadPair * Pairs;
		Library * Libs;
		ContigNode * contigs;

		int numPairs;
		int numLibs;
		int numContigs;

		int maxOverlap;
		int maxDegree;
		int minContig;
		int maxK;
};

#endif

