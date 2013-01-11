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

#ifndef DIR_NODE_H
#define DIR_NODE_H

#include <fstream>

struct DirEdge
{
	int to;
	int from;
	int flow;
};

struct DirWeight
{
	int name;
	int weight;
};

class DirNode
{
	public:

		DirNode();
		~DirNode();

		void init(int);
		void add_edge(DirEdge &);
		void fill_que(int, int *, DirNode *, int &);
		void fill_reachable(int, DirWeight *, DirNode *, int &);
		void write2file(std::ofstream & fh);
		void reset();
		void reset_flow();

		DirEdge ** edges;
		int numFilled;

		int parent;
		int num;
		bool enqued;
		bool forbidden;
		int weight;
};

#endif
