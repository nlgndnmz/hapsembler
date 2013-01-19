/*****************************************************************************
    $Author: Nilgun Donmez $
    $Date: 2012-08-25 17:04:15 -0400 (Sat, 25 Aug 2012) $

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

#ifndef CONTIG_NODE_H
#define CONTIG_NODE_H

#include <fstream>

#include "BiNode.h"
#include "ContigEdge.h"

class Containee
{
	public:
		Containee(){};
		Containee(Containee * ce)
		{
			this->contig = ce->contig;
			this->mean = ce->mean;
			this->variance = ce->variance;
			this->support = ce->support;
		}
		~Containee(){};

		int contig;
		int mean;
		int variance;
		int support;
};

class ContigNode: public BiNode<ContigEdge>
{
	public:
		ContigNode();
		~ContigNode();
		void fill_que(ContigNode *, int *, int &, bool *);
		void prune_edges(ContigNode *, int);
		void prune_nodes(ContigNode *, int);
		void remove(ContigNode *);

		ContigEdge * search(Direction, Direction, int);
		bool add_arc(ContigEdge *, int);
		bool add_fake_arc(ContigEdge *, int);

		void finalize(int, ContigNode *, Library *);
		void print_edges(std::ofstream &);
		void remove_conflicts(int);

		double coverage;
		bool forward;
		int component;

		int parent;
		bool same;
		int offset;
		int containedLen;

		LinkedList<Containee> contains;
};

ContigNode::ContigNode()
{
	coverage = 0.0;
	forward = true;
	parent = 0;
	same = true;
	offset = 0;
	containedLen = 0;
	component = 0;
}

ContigNode::~ContigNode()
{
	contains.set_delete();
	contains.clr_list();
	remove_all();	// base class method
}

void ContigNode::fill_que(ContigNode * contigs, int * que, int & end, bool * assigned)
{
	LinkedIter<ContigEdge> iter(outArcs);
	for(int d=0; d<2; d++)
	{
		for(ContigEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
		{
			int i = iter.get_flg();
			int B = e_ptr->toRead[i];

			if(!assigned[B])
			{
				que[end++] = B;
				assigned[B] = true;

				if( (d==1 && (e_ptr->arrow[i] & OUTIE)) || (d==0 && (e_ptr->arrow[i] & INNIE)) )
					contigs[B].forward = this->forward;
				else
					contigs[B].forward = !(this->forward);
			}
		}
		iter.change_list(inArcs);
	}
}

void ContigNode::prune_edges(ContigNode * contigs, int minSupport)
{
	LinkedIter<ContigEdge> iter(outArcs);
	for(int d=0; d<2; d++)
	{
		if(d) { iter.change_list(inArcs); }

		for(ContigEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
		{
			if(e_ptr->support < minSupport)
				e_ptr->flag(REMOVABLE);
		}
	}
}

void ContigNode::prune_nodes(ContigNode * contigs, int maxDeg)
{
	if(get_degree() == maxDeg)
	{
		LinkedIter<ContigEdge> iter(outArcs);
		for(int d=0; d<2; d++)
		{
			if(d) { iter.change_list(inArcs); }
			for(ContigEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
				e_ptr->flag(REMOVABLE);
		}
	}
}

void ContigNode::remove(ContigNode * contigs)
{
	LinkedIter<ContigEdge> iter(outArcs);
	for(int d=0; d<2; d++)
	{
		if(d) { iter.change_list(inArcs); }

		for(ContigEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
		{
			int i = iter.get_flg();
			int B = e_ptr->toRead[i];

			contigs[B].remove_neighbour(e_ptr->toRead[(i+1)%2]);
		}
	}
	clear_edges();
}

ContigEdge * ContigNode::search(Direction dirA, Direction dirB, int B)
{
	LinkedIter<ContigEdge> iter(get_edge_list(dirA));

	for(ContigEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
	{
		int i = iter.get_flg();
		int C = e_ptr->toRead[i];

		if(C == B && e_ptr->arrow[i] & dirB)
			return e_ptr;
	}
	return NULL;
}

bool ContigNode::add_arc(ContigEdge * oe, int which)
{
	Direction dir = (oe->arrow[(which+1)%2] & OUTIE) ? INNIE : OUTIE;
	ContigEdge * ce = search(dir, oe->arrow[which], oe->toRead[which]);

	if(ce != NULL)
	{
		ce->merge(oe);
		return false;
	}

	if(oe->arrow[(which+1)%2] & INNIE)
		inArcs.add_entry(oe, which);
	else
		outArcs.add_entry(oe, which);

	return true;
}

bool ContigNode::add_fake_arc(ContigEdge * oe, int which)
{
	Direction dir = (oe->arrow[(which+1)%2] & OUTIE) ? INNIE : OUTIE;
	ContigEdge * ce = search(dir, oe->arrow[which], oe->toRead[which]);

	if(ce != NULL)
	{
		ce->arrow[0] -= REMOVABLE;		// unflag the edge
		ce->arrow[1] -= REMOVABLE;
		return false;
	}

	if(oe->arrow[(which+1)%2] & INNIE)
		inArcs.add_entry(oe, which);
	else
		outArcs.add_entry(oe, which);

	return true;
}

void ContigNode::remove_conflicts(int maxDegree)
{
	if(this->get_degree() > maxDegree)
		this->flag |= REMOVABLE;

	LinkedIter<ContigEdge> iter(outArcs);
	for(int d=0; d<2; d++)
	{
		if(d) { iter.change_list(inArcs); }

		for(ContigEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
		{
			int f = iter.get_flg();
			int B = e_ptr->toRead[f];

			if(this->flag & REMOVABLE)	// then remove all edges
				e_ptr->flag(REMOVABLE);
			else if(f==1)				// to avoid redundancy
			{
				ContigEdge * ce[4];
				int votes[4];
				ce[0] = search(OUTIE, OUTIE, B);
				ce[1] = search(OUTIE, INNIE, B);
				ce[2] = search(INNIE, OUTIE, B);
				ce[3] = search(INNIE, INNIE, B);

				int maxSupport = 0;
				for(int y=0; y<4; y++)
				{
					if(ce[y] != NULL)
					{
						votes[y] = ce[y]->support;
						if(votes[y] > maxSupport)
							maxSupport = votes[y];
					}
				}

				int winner = 0;
				for(int y=0; y<4; y++)
				{
					if(ce[y] != NULL)
					{
						votes[y] = ce[y]->support;
						if(votes[y] < maxSupport)
							ce[y]->flag(REMOVABLE);
						else
							winner++;
					}
				}

				if(winner > 1)	// if there is a tie kill'em all
				{
					for(int y=0; y<4; y++)
					{
						if(ce[y] != NULL)
							ce[y]->flag(REMOVABLE);
					}
				}
			}
		}
	}
}

void ContigNode::finalize(int minSupport, ContigNode * contigs, Library * libs)
{
	LinkedIter<ContigEdge> iter(outArcs);
	for(int d=0; d<2; d++)
	{
		if(d) { iter.change_list(inArcs); }
		for(ContigEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
		{
			if(iter.get_flg() == 1)
				e_ptr->finalize(minSupport);
		}
	}
}

void ContigNode::print_edges(std::ofstream & fh)
{
	LinkedIter<ContigEdge> iter(outArcs);

	for(int d=0; d<2; d++)
	{
		for(ContigEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
		{
			if(iter.get_flg() == 1)
				e_ptr->print_edge(fh);
		}
		iter.change_list(inArcs);
	}
}

#endif // CONTIG_NODE_H
