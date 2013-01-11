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

#ifndef EDGE_TUPLE_H
#define EDGE_TUPLE_H

#include <fstream>

#include "GraphDef.h"

// a helper class for PathEdge
class EdgeTuple
{
	public:

		void * operator new(size_t);
		void operator delete(void *);

		EdgeTuple();
		EdgeTuple(const EdgeTuple &);
		EdgeTuple(int, Direction, EdgeLength);
		~EdgeTuple();

		void pickle(std::ofstream &);
		void unpickle(std::ifstream &);

		int read;
		Direction orientation;
		EdgeLength offset;

		static void purge();

	private:

		static void init();

		static EdgeTuple ** pool;

		static int poolnum;
		static int chunkctr;

		const static int chunksize = 1048576;	// 1024x1024
		const static int poolsize = 1048576;	// altogether this will entail a terrabyte of items
};

#endif
