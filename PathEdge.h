/*****************************************************************************
    $Author: Nilgun Donmez $
    $Date: 2011-06-10 13:52:44 -0400 (Fri, 10 Jun 2011) $

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

#ifndef PATH_EDGE_H
#define PATH_EDGE_H

#include <fstream>

#include "LinkedIter.h"
#include "EdgeTuple.h"
#include "GraphDef.h"
#include "BiEdge.h"

class PathEdge : public BiEdge
{
	public:

		void * operator new(size_t);
		void operator delete(void *);

        PathEdge();
        ~PathEdge();
        PathEdge(int, int, ArrowType, ArrowType, PathLength);

        void print_edge(std::ofstream &, int);
        bool print_edge(std::ofstream &, int, int *, int);
        int path_size();
        void add_path(LinkedList<EdgeTuple> &, int);
        bool break_tie();

        void copy_ex_first(LinkedList<EdgeTuple> &, int);
        void copy_ex_last(LinkedList<EdgeTuple> &, int);

        void pickle(std::ofstream &);
        void unpickle(std::ifstream &);

        LinkedList<EdgeTuple> path[2];
        PathLength len;

        const static ClassID data_id = PATH_EDGE;

		static void purge();

	private:

		static void init();

		static PathEdge ** pool;

		static int poolnum;
		static int chunkctr;
};

#endif // PATH_EDGE_H

