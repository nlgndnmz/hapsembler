/*****************************************************************************

	Part of Hapsembler package. See the README file for more information.
    Copyright (C) 2011-2013,  Nilgun Donmez <nild@cs.toronto.edu>

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

#ifndef READ_NODE_H
#define READ_NODE_H

#include <fstream>

#include "Containee.h"
#include "GraphDef.h"
#include "BiNode.h"
#include "OvlEdge.h"

// the main component of the OvlGraph
class ReadNode: public BiNode<OvlEdge>
{
    public:
        ReadNode();
        ~ReadNode();

        OvlEdge * search_node(int, Direction, Direction);

        void pickle(std::ofstream &);
        void unpickle(std::ifstream &, int, int, ReadNode *);

        void process(int, SDNode *, int *, int, int &, Direction, EdgeLength);

        // --- start of Containee related ---
        void print_contains(std::ofstream &);
        void add_containee(int, bool, EdgeLength, EdgeLength);
        int get_parents(int *);

        bool isChild();
        bool isParent();
        int num_kids();
        void reset_parent();

        int parent;
        bool same;
        EdgeLength offset;
        LinkedList<Containee> contains;
        // --- end of Containee related ---

        const static ClassID data_id = READ_NODE;		// like version id, but for data structures instead
};

ReadNode::ReadNode()
{
	parent = 0;
	same = true;
	offset = 0;
	contains.set_delete();
}

ReadNode::~ReadNode()
{
	// nothing to delete
}

// pickle the ReadNode to the file given
void ReadNode::pickle(std::ofstream & fh)
{
	int num = get_degree();
	fh << (int)data_id <<" "<< num <<" "<< id <<" "<< len <<" "<< flag <<" "
		<< parent << " " << (int)same <<" "<< offset <<" "<< contains.get_size() << "\n";

	LinkedIter<Containee> iter(contains);
	for(Containee * ct = iter.get_first(); ct != NULL; ct = iter.get_next())
		ct->pickle(fh);

	BiNode<OvlEdge>::pickle(fh);
}

// unpickle the node from file
void ReadNode::unpickle(std::ifstream & fh, int num, int id, ReadNode * nodes)
{
	int s, o, n;
	this->id = id;
	fh >> len >> flag >> parent >> s >> o >> n;
	same = (s==0) ? false : true;
	offset = (EdgeLength) o;

	for(int i=0; i<n; i++)
	{
		Containee * ct = new Containee();
		ct->unpickle(fh);
		contains.add_entry(ct);
	}

	for(int i=0; i<num; i++)
	{
        OvlEdge * oe = new OvlEdge();
		oe->unpickle(fh);
		int k = (oe->toRead[0] == id) ? 1 : 0;		// nodes should not have self-loops
		if(id < oe->toRead[k] || (id == oe->toRead[k] && k==1))
		{
            add_arc(oe, k);
            nodes[oe->toRead[k]].add_arc(oe, rev(k));
		}
		else
            delete oe;
	}
}

// used by verify_distance (added so that contained reads can be processed as well)
void ReadNode::process(int A, SDNode * seen, int * seenlist, int target, int & counter, Direction dir, EdgeLength val)
{
	if(seen[A].target != target)
	{
		seen[A].target = target;
		seen[A].in = -1;
		seen[A].out = -1;
		seenlist[counter++] = A;

		LinkedIter<Containee> iter(contains);
		for(Containee * ct = iter.get_first(); ct != NULL; ct = iter.get_next())
		{
			seen[ct->read].target = target;
			seen[ct->read].in = -1;
			seen[ct->read].out = -1;
			seenlist[counter++] = ct->read;
		}
	}

	if(dir == OUTIE)
		seen[A].out = val;
	else
		seen[A].in = val;

	LinkedIter<Containee> iter(contains);
	for(Containee * ct = iter.get_first(); ct != NULL; ct = iter.get_next())
	{
		if(dir == OUTIE)
		{
			if(ct->same)
				seen[ct->read].out = val + ct->end;
			else
				seen[ct->read].in = val + ct->begin;
		}
		else
		{
			if(ct->same)
				seen[ct->read].in = val + ct->begin;
			else
				seen[ct->read].out = val + ct->end;
		}
	}
}

// get the list of nodes containing this node (i.e. nodes with non-proper overlaps)
int ReadNode::get_parents(int * ptr)
{
	int counter = 0;

	LinkedIter<OvlEdge> iter(inArcs);
	for(OvlEdge * oe = iter.get_first(); oe != NULL; oe = iter.get_next())
	{
		int i = iter.get_flg();
		if(oe->arrow[i] & NON_PROPER)	// a containing read
			ptr[counter++] = oe->toRead[i];
	}
	return counter;
}

// add a new read contained by this read
void ReadNode::add_containee(int read, bool same, EdgeLength begin, EdgeLength end)
{
	Containee * ct = new Containee(read, same, begin, end);
	contains.add_entry(ct);
}

bool ReadNode::isChild()
{
	return (parent != 0);
}

bool ReadNode::isParent()
{
	return (contains.get_size() > 0);
}

int ReadNode::num_kids()
{
	return contains.get_size();
}

void ReadNode::reset_parent()
{
	parent = 0;
	same = true;
	offset = 0;
}

void ReadNode::print_contains(std::ofstream & fh)
{
	LinkedIter<Containee> iter(contains);
	for(Containee * ct = iter.get_first(); ct != NULL; ct = iter.get_next())
		fh << "-2 " << ct->read <<" "<< (int)(ct->same) <<" "<< (int)(ct->begin) <<" "<< (int)(ct->end) << "\n";
}

#endif // READ_NODE_H

