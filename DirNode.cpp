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

#include <fstream>

#include "DirNode.h"

DirNode::DirNode()
{
	edges = 0;
	numFilled = 0;
	parent = -1;
	num = -1;
	enqued = false;
	forbidden = false;
	weight = 0;
}

void DirNode::init(int num)
{
	edges = new DirEdge*[num];
	numFilled = 0;
}

void DirNode::reset()
{
	numFilled = 0;
	parent = -1;
	num = -1;
	enqued = false;
	forbidden = false;
}

DirNode::~DirNode()
{
	delete [] edges;
}

void DirNode::add_edge(DirEdge & de)
{
	edges[numFilled++] = (&de);
}

void DirNode::fill_que(int ID, int * que, DirNode * dirnodes, int & count)
{
	int B = -1;
	for(int i=0; i<numFilled; i++)
	{
		bool ok = false;
		if(edges[i]->from == ID && edges[i]->flow == 0)		// then the residual is (capacity - flow)
		{
			B = edges[i]->to;
			ok = true;
		}
		else if(edges[i]->to == ID && edges[i]->flow == 1)		// means we can push back flow
		{
			B = edges[i]->from;
			ok = true;
		}

		if(ok && !dirnodes[B].enqued && !dirnodes[B].forbidden)
		{
			que[count++] = B;
			dirnodes[B].parent = ID;
			dirnodes[B].num = i;
			dirnodes[B].enqued = true;
		}
	}
}

void DirNode::fill_reachable(int ID, DirWeight * que, DirNode * dirnodes, int & count)
{
	for(int i=0; i<numFilled; i++)
	{
		if(edges[i]->from == ID && edges[i]->flow == 1)
		{
			int B = edges[i]->to;
			if(!dirnodes[B].enqued)
			{
				dirnodes[B].enqued = true;
				que[count].name = B;
				que[count++].weight = dirnodes[B].weight;
			}
		}
	}
}

void DirNode::reset_flow()
{
	for(int i=0; i<numFilled; i++)
		edges[i]->flow = 0;
}

void DirNode::write2file(std::ofstream & fh)
{
	for(int i=0; i<numFilled; i++)
		fh << "from: " << edges[i]->from << " to: " << edges[i]->to << " flow: " << edges[i]->flow << "\n";
}
