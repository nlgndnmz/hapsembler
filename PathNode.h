/*****************************************************************************
    $Author: Nilgun Donmez $
    $Date: 2012-02-02 16:56:38 -0500 (Thu, 02 Feb 2012) $

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

#ifndef PATH_NODE_H
#define PATH_NODE_H

#include <fstream>

#include "BiNode.h"
#include "PathEdge.h"

// this class can be used both for MatePair and ReadNode related paths
class PathNode: public BiNode<PathEdge>
{
    public:

        PathNode() { comp = 0; count = 0; forward = false; };
        ~PathNode() {};

        void print_edges(std::ofstream &, int *, int);
        void print_edges(std::ofstream &, PathNode *);
        void print_edges_ss(std::ofstream &, int *, int);
        void pickle(std::ofstream &);
        void unpickle(std::ifstream &, int, int, PathNode *);

        const static ClassID data_id = PATH_NODE;

        int comp;
        bool forward;
        int count;
};

// if reporting for contigs after mate dfs, report paths with at least 3 unvisited nodes, otherwise report everything
void PathNode::print_edges(std::ofstream & fh, int * reported = NULL, int minimum = 3)
{
	LinkedIter<PathEdge> iter(outArcs);
	for(int d=0; d<2; d++)
	{
		if(d)	iter.change_list(inArcs);
		for(PathEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
		{
			int i = iter.get_flg();
			if(id < e_ptr->toRead[i] || (id == e_ptr->toRead[i] && i==1))	// to write the edge only once
			{
				if(reported == NULL)
				{
					e_ptr->print_edge(fh, i);
					fh << "-1 0" << "\n";
				}
				else if( e_ptr->print_edge(fh, i, reported, minimum) )
					fh << "-1 0" << "\n";
			}
		}
	}
}

void PathNode::print_edges(std::ofstream & fh, PathNode * paths)
{
	LinkedIter<PathEdge> iter(outArcs);
	for(int d=0; d<2; d++)
	{
		if(d)	iter.change_list(inArcs);
		for(PathEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
		{
			int i = iter.get_flg();
			if(id < e_ptr->toRead[i] || (id == e_ptr->toRead[i] && i==1))	// to write the edge only once
			{
				if(d==1 && this->get_in_degree() == 1)
					fh << "0 " << OUTIE << " " << this->id << " " << this->len << "\n";
				else if(d==0 && this->get_out_degree() == 1)
					fh << "0 " << INNIE << " " << this->id << " " << this->len << "\n";

				LinkedIter<EdgeTuple> iter(e_ptr->path[i]);
				EdgeTuple * et = iter.get_first();
				EdgeTuple * et2 = iter.get_next();
				while(et2 != NULL)
				{
					fh << et->offset << " " << (et->orientation & (INNIE|OUTIE)) << " " << et->read << " " << paths[et->read].len << "\n";
					et = et2;
					et2 = iter.get_next();
				}

				if((e_ptr->arrow[i] & INNIE) && paths[et->read].get_in_degree() == 1)
					fh << et->offset << " " << INNIE << " " << et->read << " " << paths[et->read].len << "\n";
				else if((e_ptr->arrow[i] & OUTIE) && paths[et->read].get_out_degree() == 1)
					fh << et->offset << " " << OUTIE << " " << et->read << " " << paths[et->read].len << "\n";

				fh << "-1 0" << "\n";
			}
		}
	}
}

// for strand specific reporting
void PathNode::print_edges_ss(std::ofstream & fh, int * reported = NULL, int minimum = 3)
{
	LinkedIter<PathEdge> iter(outArcs);
	for(PathEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
	{
		int i = iter.get_flg();
		if(reported == NULL)
		{
			e_ptr->print_edge(fh, i);
			fh << "-1 0" << "\n";
		}
		else if( e_ptr->print_edge(fh, i, reported, minimum) )
			fh << "-1 0" << "\n";
	}
}

// pickle the PathNode to the file given
void PathNode::pickle(std::ofstream & fh)
{
	int num = get_degree();
	fh << (int)data_id <<" "<< num <<" "<< id <<" "<< len <<" "<< (int) flag << "\n";
	BiNode<PathEdge>::pickle(fh);	// call the base class method
}

// unpickle the PathNode from file
void PathNode::unpickle(std::ifstream & fh, int num, int id, PathNode * nodes)
{
	this->id = id;
	fh >> len >> flag;

	for(int i=0; i<num; i++)
	{
        PathEdge * pe = new PathEdge();
		pe->unpickle(fh);
		int k = (pe->toRead[0] == id) ? 1 : 0;		// nodes should not have self-loops
		if(id < pe->toRead[k] || (id == pe->toRead[k] && k==1))
		{
            add_arc(pe, k);
            nodes[pe->toRead[k]].add_arc(pe, rev(k));
		}
		else
            delete pe;
	}
}

#endif // PATH_NODE_H
