/*****************************************************************************
    $Author: $
    $Date: $

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

#ifndef U_GRAPH_H
#define U_GRAPH_H

#include "DirNode.h"

class UNode;

class UGraph
{
	public:
		UGraph(int, int, int);
		~UGraph();

		void add_edge(int, int);
		void add_node(int, int);
		int oddcycle(int *, int, bool &);

	private:

		void add_diredge(int, int);
		void augment_X_edges(int *, int);
		void augment_Y_edges(int *, int);
		void augment_S_edges(int *, int, int *, int);

		int ford_fulkerson(int, int *);
		int solve_cut(int, int *, int *, DirWeight *, int);
		void clear_dirnodes(int);
		void repartition(int *, int &, int *, int &, int *);

		int numNodes;
		int numEdges;

		UNode * nodes;

		DirNode * dirnodes;
		DirEdge * diredges;

};

#endif
