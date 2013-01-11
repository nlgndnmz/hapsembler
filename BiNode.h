/*****************************************************************************
    $Author: Nilgun Donmez $
    $Date: 2012-01-18 20:48:51 -0500 (Wed, 18 Jan 2012) $

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

#ifndef BI_NODE_H
#define BI_NODE_H

#include <fstream>

#include "LinkedIter.h"
#include "GraphDef.h"

// a templated class for a node in a bidirected graph
template <class EDGE_TYPE>
class BiNode
{
    public:

		// methods

        BiNode();
        ~BiNode();

        int get_in_degree();
        int get_out_degree();
        int get_degree();

        void add_arc(EDGE_TYPE *, int);

     	void clear_edges();
     	void sort_edges();
     	bool is_end_node();
     	bool is_single_end();
     	bool is_middle();

     	bool isNeighbour(int, Direction, Direction);

     	int remove_marked(EdgeFlag);
     	int remove_marked(BiNode<EDGE_TYPE> *, NodeFlag);
     	int remove_unmarked(EdgeFlag);
     	int remove_neighbour(int);
     	int remove_all();

     	int get_neighbours(Direction, int *, bool);
     	LinkedList<EDGE_TYPE> & get_edge_list(Direction);

     	void pickle(std::ofstream &);

     	// data members

		int id;
		NodeLength len;
		NodeFlag flag;

        LinkedList<EDGE_TYPE> inArcs;
        LinkedList<EDGE_TYPE> outArcs;
};

// initialzes the primitive types and also allocates the edge lists
template <class EDGE_TYPE>
BiNode<EDGE_TYPE>::BiNode()
{
	id = 0;
	len = 0;
	flag = NONE;

	inArcs.set_no_delete();
	outArcs.set_no_delete();
}

// first deletes the contents of the edge lists
template <class EDGE_TYPE>
BiNode<EDGE_TYPE>::~BiNode()
{
	clear_edges();
}

// returns the indegree
template <class EDGE_TYPE>
int BiNode<EDGE_TYPE>::get_in_degree()
{
	return inArcs.get_size();
}

// returns the outdegree
template <class EDGE_TYPE>
int BiNode<EDGE_TYPE>::get_out_degree()
{
	return outArcs.get_size();
}

// returns the total degree of the node
template <class EDGE_TYPE>
int BiNode<EDGE_TYPE>::get_degree()
{
	return ( inArcs.get_size() + outArcs.get_size() );
}

// deletes all the elements in the edge lists and flags the node as dead (does not actually delete lists themselves)
template <class EDGE_TYPE>
void BiNode<EDGE_TYPE>::clear_edges()
{
	flag |= DEAD;

	inArcs.set_delete();
	inArcs.clr_list();

	outArcs.set_delete();
	outArcs.clr_list();
}

// gets the edge list in opposite direction to orient
template <class EDGE_TYPE>
LinkedList<EDGE_TYPE> & BiNode<EDGE_TYPE>::get_edge_list(Direction orient)
{
	if(orient & INNIE)
		return outArcs;
	return inArcs;
}

// adds an arc to the node
template <class EDGE_TYPE>
void BiNode<EDGE_TYPE>::add_arc(EDGE_TYPE * oe, int which)
{
	if(oe->arrow[(which+1)%2] & INNIE)
		inArcs.add_entry(oe, which);
	else
		outArcs.add_entry(oe, which);
}

// sorts the edges based on the length of the edges
template <class EDGE_TYPE>
void BiNode<EDGE_TYPE>::sort_edges()
{
	inArcs.sort();
	outArcs.sort();
}

// removes edges marked with the flag F
template <class EDGE_TYPE>
int BiNode<EDGE_TYPE>::remove_marked(EdgeFlag F)
{
	int counter = 0;
	LinkedIter<EDGE_TYPE> iter(inArcs);
	for(int d=0; d<2; d++)
	{
		if(d)	iter.change_list(outArcs);

		EDGE_TYPE * oe = iter.get_first();
		while(oe != NULL)
		{
			if(oe->flagged(F))		// we have a hit
			{
				counter += 1;
				if(oe->flagged(OK))			// if the edge has been deleted by the other node as well
					oe = iter.del_cur(true);	// does real deletion
				else
				{
					oe->flag(OK);
					oe = iter.del_cur();	// delete the LinkedNode but do not delete the actual data
				}
			}
			else
			{
				oe = iter.get_next();	// this will account for the fact that we deleted the current one
			}
		}
	}
	return counter;
}

// removes neighbouring nodes marked with the flag F
template <class EDGE_TYPE>
int BiNode<EDGE_TYPE>::remove_marked(BiNode<EDGE_TYPE> * nd_ptr, NodeFlag F)
{
	int counter = 0;
	LinkedIter<EDGE_TYPE> iter(inArcs);
	for(int d=0; d<2; d++)
	{
		if(d)	iter.change_list(outArcs);

		EDGE_TYPE * oe = iter.get_first();
		while(oe != NULL)
		{
			if(nd_ptr[ oe->toRead[iter.get_flg()] ].flag & F)		// we have a hit
			{
				counter += 1;
				if(oe->flagged(OK))			// if the edge has been deleted by the other node as well
					oe = iter.del_cur(true);	// does real deletion
				else
				{
					oe->flag(OK);
					oe = iter.del_cur();	// delete the LinkedNode but do not delete the actual data
				}
			}
			else
			{
				oe = iter.get_next();	// this will account for the fact that we deleted the current one
			}
		}
	}
	return counter;
}

// removes edges if they are not marked with the flag F
template <class EDGE_TYPE>
int BiNode<EDGE_TYPE>::remove_unmarked(EdgeFlag F)
{
	int counter = 0;
	LinkedIter<EDGE_TYPE> iter(inArcs);
	for(int d=0; d<2; d++)
	{
		if(d)	iter.change_list(outArcs);

		EDGE_TYPE * oe = iter.get_first();
		while(oe != NULL)
		{
			if(!(oe->flag & F))		// we have a hit
			{
				counter += 1;
				if(oe->flagged(OK))			// if the edge has been deleted by the other node as well
					oe = iter.del_cur(true);	// does real deletion
				else
				{
					oe->flag(OK);
					oe = iter.del_cur();	// delete the LinkedNode but do not delete the actual data
				}
			}
			else
			{
				oe = iter.get_next();	// this will account for the fact that we deleted the current one
			}
		}
	}
	return counter;
}

// remove all edges
template <class EDGE_TYPE>
int BiNode<EDGE_TYPE>::remove_all()
{
	int counter = 0;
	LinkedIter<EDGE_TYPE> iter(inArcs);
	for(int d=0; d<2; d++)
	{
		if(d)	iter.change_list(outArcs);

		EDGE_TYPE * oe = iter.get_first();
		while(oe != NULL)
		{
			counter += 1;
			if(oe->flagged(OK))			// if the edge has been deleted by the other node as well
				oe = iter.del_cur(true);	// does real deletion
			else
			{
				oe->flag(OK);
				oe = iter.del_cur();	// delete the LinkedNode but do not delete the actual data
			}
		}
	}
	return counter;
}

// removes the neighbour with id u
template <class EDGE_TYPE>
int BiNode<EDGE_TYPE>::remove_neighbour(int u)
{
	int counter = 0;
	LinkedIter<EDGE_TYPE> iter(inArcs);
	for(int d=0; d<2; d++)
	{
		if(d)	iter.change_list(outArcs);

		EDGE_TYPE * oe = iter.get_first();
		while(oe != NULL)
		{
			int i = iter.get_flg();
			if(oe->toRead[i] == u)		// we have a hit!
			{
				counter += 1;
                if(oe->flagged(OK))
                    oe = iter.del_cur(true);
                else
                {
                    oe->flag(OK);
                    oe = iter.del_cur();
                }
			}
			else
			{
				oe = iter.get_next();
			}
		}
	}
	return counter;
}

// removes the neighbour with id u
template <class EDGE_TYPE>
bool BiNode<EDGE_TYPE>::isNeighbour(int u, Direction dir1, Direction dir2)
{
	LinkedIter<EDGE_TYPE> iter(inArcs);
	if(dir1 & OUTIE)
		iter.change_list(outArcs);

	for(EDGE_TYPE * oe = iter.get_first(); oe != NULL; oe = iter.get_next())
	{
		int i = iter.get_flg();
		if(oe->toRead[i] == u && oe->arrow[i] & dir2)		// we have a hit!
			return true;
	}
	return false;
}

// returns true if the node has only incoming or only outgoing edges
template <class EDGE_TYPE>
bool BiNode<EDGE_TYPE>::is_end_node()
{
	int a = inArcs.get_size();
	int b = outArcs.get_size();

	if(a+b > 0 && (a == 0 || b == 0))
		return true;

	return false;
}

// returns true if the node has only one edge
template <class EDGE_TYPE>
bool BiNode<EDGE_TYPE>::is_single_end()
{
	if(get_degree() == 1)
		return true;

	return false;
}

// returns true if the node has exactly one incoming and one outgoing edge
template <class EDGE_TYPE>
bool BiNode<EDGE_TYPE>::is_middle()
{
	if(get_in_degree() == 1 && get_out_degree() == 1)
		return true;
	return false;
}

// gets a list of all the neighbours, if noparent is true do not get neighbours with non-proper edges
template <class EDGE_TYPE>
int BiNode<EDGE_TYPE>::get_neighbours(Direction in_out, int * ptr, bool noparent)
{
	int counter = 0;
	int i;

	if(in_out & INNIE)
	{
		LinkedIter<EDGE_TYPE> iter(inArcs);
		for(EDGE_TYPE * oe = iter.get_first(); oe != NULL; oe = iter.get_next())
		{
			i = iter.get_flg();
			if(noparent && (oe->arrow[i] & NON_PROPER))	;
			else
				ptr[counter++] = (oe->arrow[i] & INNIE) ? oe->toRead[i] : -(oe->toRead[i]);
		}
	}

	if(in_out & OUTIE)
	{
		LinkedIter<EDGE_TYPE> iter(outArcs);
		for(EDGE_TYPE * oe = iter.get_first(); oe != NULL; oe = iter.get_next())
		{
			i = iter.get_flg();
			if(noparent && (oe->arrow[i] & NON_PROPER))	;
			else
				ptr[counter++] = (oe->arrow[i] & INNIE) ? oe->toRead[i] : -(oe->toRead[i]);
		}
	}
	return counter;
}

// pickles the edges of the node
template <class EDGE_TYPE>
void BiNode<EDGE_TYPE>::pickle(std::ofstream & fh)
{
	LinkedIter<EDGE_TYPE> iter(outArcs);
	for(int d=0; d<2; d++)
	{
		if(d)	iter.change_list(inArcs);
		for(EDGE_TYPE * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
            e_ptr->pickle(fh);
	}
}

#endif // BI_NODE_H

