/*****************************************************************************
    $Author: Nilgun Donmez $
    $Date: 2011-12-08 18:56:55 -0500 (Thu, 08 Dec 2011) $

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

#ifndef OVL_EDGE_H
#define OVL_EDGE_H

#include <fstream>

#include "GraphDef.h"
#include "BiEdge.h"

// a class for bidirected overlap edges
class OvlEdge: public BiEdge
{
	public:

		void * operator new(size_t);
		void operator delete(void *);

		OvlEdge();
		OvlEdge(int, int, EdgeLength, EdgeLength, ArrowType, ArrowType);
		~OvlEdge();

		void print_edge(std::ofstream &, int);
		void pickle(std::ofstream &);
		void unpickle(std::ifstream &);

		EdgeLength len[2];

		static void shrink();
		static void purge();

	private:

		bool alive;

		static void init();

		static OvlEdge ** pool;

		static int poolnum;
		static int chunkctr;
};

#endif // OVL_EDGE_H
