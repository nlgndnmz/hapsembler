/*****************************************************************************
    $Author: Nilgun Donmez $
    $Date: 2012-02-15 22:44:52 -0500 (Wed, 15 Feb 2012) $

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

#ifndef CONTIG_EDGE_H
#define CONTIG_EDGE_H

#include <fstream>

#include "GraphDef.h"
#include "BiEdge.h"

class ContigEdge: public BiEdge
{
	public:

		void * operator new(size_t);
		void operator delete(void *);

		ContigEdge(int, int, ArrowType, ArrowType, int, int, int, double);
		explicit ContigEdge(const ContigEdge &);
		ContigEdge();
		~ContigEdge();

		void merge(ContigEdge *);
		void finalize(int);
		void print_edge(std::ofstream &);

		int support;
		int len[2];
		int lib;
		double weight;

		static void purge();
		static void shrink();

	private:

		static void init();

		static ContigEdge ** pool;

		static int poolnum;
		static int chunkctr;
};

#endif // CONTIG_EDGE_H
