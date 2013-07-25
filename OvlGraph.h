/*****************************************************************************
    $Author: Nilgun Donmez $
    $Date: 2012-08-28 18:36:05 -0400 (Tue, 28 Aug 2012) $

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

#ifndef OVL_GRAPH_H
#define OVL_GRAPH_H

#include <fstream>

#ifdef USEOPENMP
	#include <omp.h>
#endif

#include "GraphDef.h"

class PathNode;
class MatePair;
class ReadNode;

class OvlGraph
{
	public:
		OvlGraph(int, bool);
		~OvlGraph();

		// populators
		void add_node(int, int);
		void add_pair(int, short int);
		void add_arc(int, int, int, int, int, int, int, bool);

		void build_pair_graph();

		// twin methods
		int remove_contained_reads_alt();
		void remove_contained_reads();
		void remove_contained_mates();

		void reduce_read_graph();
		void reduce_mate_graph();

		int remove_read_bubbles();
		int remove_mate_bubbles();

		int remove_read_buds();
		int remove_mate_buds();

		void read_chain_collapse();
		void mate_chain_collapse();

		int check_read_graph();
		int check_mate_graph();

		void sort_read_edges();
		void sort_mate_edges();

		void print_read_contigs(char *);
		int print_mate_contigs(char *);

		void simplify_read_chains(int);
		void simplify_mate_chains(int);

		// simplifiers
		long int venom();
		void bridge_matepairs();

		// settlers
		void set_pairs(int, int);
		void set_flexi(double);

		// miscellaneous methods
		int check_consistency();
		void topsort(char *);
		void remove_short_mate_overlaps(int);
		void remove_short_read_overlaps(int);
		void mark_single_nodes(double, double);

		void calibrate_libraries(char *);

		Library * Libs;

	private:

		bool checkPath(int, int &, int *);

		int verify_distance(int , int , int , QStruct * , SDNode * , PStruct * , Direction, int *, int &, SDNode *, bool);
		int close_gap(std::ofstream & , int , int , int , QStruct * , SDNode * , PStruct * , int * , int * , int, Direction, Direction, int *, bool &);

		bool isConsistent(int , int , int , QStruct * , SDNode * , PStruct * , int, Direction, Direction);

		// data members
		bool onestrand;

		int numNodes;
		int numPairs;
		int numLibs;

		ReadNode * Nodes;
		MatePair * Pairs;
		PathNode * Mates;
		PathNode * Reads;
};

#endif // OVL_GRAPH_H
