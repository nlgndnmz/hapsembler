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

#ifndef SCAFF_NODE_H
#define SCAFF_NODE_H

#include "LinkedIter.h"

class ScaffEdge;

class ScaffNode
{
	public:

		ScaffNode();
		~ScaffNode();

		int add_in_edge(int, int, int);
		int add_out_edge(int, int, int);
		void clear_edges();

		void finalize();

		bool isInternal();
		bool isEligible(ScaffNode *);
		bool isStart();
		bool isEnd();

		bool isLonger(int);

		int id;
		bool forward;
		int len;

		int parent;
		int offset;
		int offsetVar;

		int totalLen;
		int totalVar;

		LinkedList<ScaffEdge> inArcs;
		LinkedList<ScaffEdge> outArcs;
		LinkedList<ScaffEdge> contained;
};

#endif
