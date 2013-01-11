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

#ifndef SCAFF_EDGE_H
#define SCAFF_EDGE_H

#include "LinkedIter.h"

struct Fuzzy
{
	int mu;
	int var;
};

class ScaffEdge
{
	public:

		explicit ScaffEdge(int);
		explicit ScaffEdge(ScaffEdge *);
		~ScaffEdge();

		void add_fuzz(int, int);
		void finalize();

		int toContig;
		int mean;
		int variance;
		LinkedList<Fuzzy> fuzzies;
};

#endif
