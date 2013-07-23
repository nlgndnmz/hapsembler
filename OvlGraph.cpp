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

#include <cmath>
#include <fstream>
#include <cstring>
#include <iostream>

#include "OvlGraphTempl.h"
#include "OvlGraph.h"

using namespace std;

// Create a graph with the given number of nodes
OvlGraph::OvlGraph(int num, bool ss)
{
	numNodes = num;
	numPairs = 0;

	Nodes = new ReadNode[numNodes + 1];

	Reads = NULL;
	Pairs = NULL;
	Mates = NULL;
	Libs = NULL;

	onestrand = ss;
}

// delete the graph
OvlGraph::~OvlGraph()
{
	purge_edges<ReadNode, OvlEdge>(Nodes, numNodes);
	purge_edges<PathNode, PathEdge>(Reads, numNodes);

	purge_edges<MatePair, MateEdge>(Pairs, numPairs);
	purge_edges<PathNode, PathEdge>(Mates, numPairs);

	delete [] Nodes;
	delete [] Reads;

	LinkedNode<OvlEdge>::purge();
	LinkedNode<MateEdge>::purge();
	LinkedNode<PathEdge>::purge();
	LinkedNode<EdgeTuple>::purge();
	LinkedNode<Containee>::purge();

	if(Pairs != NULL)
		delete [] Pairs;
	if(Mates != NULL)
		delete [] Mates;
	if(Libs != NULL)
		delete [] Libs;

	OvlEdge::purge();
	MateEdge::purge();	// if no MateEdge is allocated, this will have no effect
	PathEdge::purge();
	EdgeTuple::purge();
}

void OvlGraph::remove_contained_reads()
{
	remove_contained_nodes<ReadNode, OvlEdge>(Nodes, numNodes);
}

void OvlGraph::remove_contained_mates()
{
	remove_contained_nodes<MatePair, MateEdge>(Pairs, numPairs);
}

void OvlGraph::reduce_read_graph()
{
	reduce_graph<ReadNode, OvlEdge>(Nodes, numNodes, 0);
}

void OvlGraph::reduce_mate_graph()
{
	reduce_graph<MatePair, MateEdge>(Pairs, numPairs, Nodes);
}

int OvlGraph::remove_read_bubbles()
{
	return remove_bubbles<ReadNode, OvlEdge>(Nodes, numNodes);
}

int OvlGraph::remove_mate_bubbles()
{
	return remove_bubbles<MatePair, MateEdge>(Pairs, numPairs);
}

int OvlGraph::remove_read_buds()
{
	return remove_buds<ReadNode, OvlEdge>(Nodes, numNodes);
}

int OvlGraph::remove_mate_buds()
{
	return remove_buds<MatePair, MateEdge>(Pairs, numPairs);
}

void OvlGraph::read_chain_collapse()
{
	Reads = new PathNode[numNodes + 1];
	for(int i=1; i<=numNodes; i++)
	{
		Reads[i].id = Nodes[i].id;
		Reads[i].len = Nodes[i].len;		// flag is initially zero
	}
	chain_collapse<ReadNode, OvlEdge>(Nodes, Reads, numNodes);
}

void OvlGraph::mate_chain_collapse()
{
	Mates = new PathNode[numPairs + 1];
	for(int i=1; i<=numPairs; i++)
		Mates[i].id = Pairs[i].id;
	chain_collapse<MatePair, MateEdge>(Pairs, Mates, numPairs);
}

int OvlGraph::plot_readovl_distribution(char * fn, int t)
{
	return plot_distribution<ReadNode>(fn, 1000, t, Nodes, numNodes);
}

int OvlGraph::plot_mateovl_distribution(char * fn, int t)
{
	return plot_distribution<MatePair>(fn, 1000, t, Pairs, numPairs);
}

int OvlGraph::check_read_graph()
{
	return check_graph<ReadNode>(Nodes, numNodes);
}

int OvlGraph::check_mate_graph()
{
	return check_graph<MatePair>(Pairs, numPairs);
}

void OvlGraph::sort_read_edges()
{
	for(int A=1; A<=numNodes; A++)
		Nodes[A].sort_edges();
}

void OvlGraph::sort_mate_edges()
{
	for(int A=1; A<=numPairs; A++)
		Pairs[A].sort_edges();
}

void OvlGraph::print_read_contigs(char * fn)
{
	print_contigs(fn, Reads, numNodes, onestrand);
}

void OvlGraph::simplify_read_chains(int maxRemoval)
{
	int x = 1;
	int popped = 0;
	cout << "buds removed: " << remove_buds(Reads, numNodes, maxRemoval) << endl;
	path_merge(Reads, numNodes);
	while(x > 0)
	{
		x = remove_bubbles(Reads, numNodes, maxRemoval);
		popped += x;
		path_merge(Reads, numNodes);
	}
	cout << "bubbles removed: " << popped << endl;
	path_sync<ReadNode, OvlEdge>(Nodes, numNodes, Reads);
	check_graph(Nodes, numNodes);
	check_graph(Reads, numNodes);
}

void OvlGraph::simplify_mate_chains(int maxRemoval)
{
	int x = 1, y = 1;
	int popped = 0;
	int peaked = 0;
	check_graph(Mates, numPairs);
	cout << "buds removed: " << remove_buds(Mates, numPairs, maxRemoval) << endl;
	path_merge(Mates, numPairs);
	while(x+y > 0)
	{
		x = remove_bubbles(Mates, numPairs, maxRemoval);
		popped += x;
		path_merge(Mates, numPairs);
		y = peekaboo(Mates, numPairs);
		peaked += y;
		path_merge(Mates, numPairs);
	}
	cout << "bubbles removed: " << popped << endl;
	cout << "cross edges removed: " << peaked << endl;
	check_graph(Mates, numPairs);
}

void OvlGraph::set_pairs(int num, int numL)
{
	numPairs = num;
	numLibs = numL;
	Pairs = new MatePair[numPairs + 1];
	Libs = new Library[numL+1];
}

void OvlGraph::set_flexi(double flexi)
{
	MatePair::FLEX = flexi;
}

void OvlGraph::topsort(char * filename)
{
	ofstream fh;
	open_n_check(fh, filename);

	int component = 1;
	int * que = new int[numPairs+1];
	int * S = new int[numPairs+1];
	int * L = new int[numPairs+1];

	for(int x=1; x<=numPairs; x++)
	{
		if(Mates[x].comp == 0 && Mates[x].get_degree() > 0)
		{
			int winner = matesort(Mates, x, ++component, que, S, L);
			Direction dir = (que[winner] & OUTIE) ? INNIE : OUTIE;
			while(S[winner] != 0)
			{
				LinkedIter<PathEdge> iter(Mates[winner].get_edge_list(dir));
				for(PathEdge * pe = iter.get_first(); pe != NULL; pe = iter.get_next())
				{
					int i = iter.get_flg();
					if(pe->toRead[i] == S[winner])
					{
						winner = pe->toRead[i];
						dir = (que[winner] & OUTIE) ? INNIE : OUTIE;
						fh << pe->len <<" "<< pe->path_size() <<" "<< winner << "\n";
					}
				}
			}
			fh << "-1 0" << "\n";
		}
	}
	check_n_close(fh);

	delete [] que;
	delete [] S;
	delete [] L;
}

// add a node to the graph
void OvlGraph::add_node(int ID, NodeLength length)
{
	if(ID > numNodes)
		return;

	if(length < 0 || ID < 0)
		throw "OvlGraph: Read length and ID must be positive";

	Nodes[ID].id = ID;
	Nodes[ID].len = length;		// flag is initially zero
}

// add a mate pair to graph. Note that this method should be called after all the ReadNodes are initialized.
void OvlGraph::add_pair(int ID, short int lib)
{
	if(ID > numPairs)
		return;

	Pairs[ID].lib = lib;
	Pairs[ID].id = ID;
}

// add an OvlEdge to the graph. For contained reads add two (non-proper) edges
void OvlGraph::add_arc(int first_read, int second_read, int first_start, int second_start,
	int first_end, int second_end, int complement, bool keep)
{
	if(first_read > numNodes || second_read > numNodes)		// reads are out of bound
		return;

	OvlEdge * oe1 = NULL;
	OvlEdge * oe2 = NULL;

	if( (first_end - first_start) >= Nodes[first_read].len) // first read is contained
    {
    	assert(first_start == 0);
    	assert(first_end == Nodes[first_read].len);
		Nodes[first_read].flag |= CONTAINED;

		if(!keep)
			return;

        if(complement)
        {
        	oe1 = new OvlEdge(first_read, second_read, Nodes[second_read].len - second_end, 0 - second_start, OUTIE, OUTIE|NON_PROPER);
			oe2 = new OvlEdge(first_read, second_read, second_start, second_end - Nodes[second_read].len, INNIE, INNIE|NON_PROPER);
        }
        else
        {
        	oe1 = new OvlEdge(first_read, second_read, Nodes[second_read].len - second_end, 0 - second_start, OUTIE, INNIE|NON_PROPER);
			oe2 = new OvlEdge(first_read, second_read, second_start, second_end - Nodes[second_read].len,  INNIE, OUTIE|NON_PROPER);
        }

		Nodes[first_read].add_arc(oe1, 1);
		Nodes[first_read].add_arc(oe2, 1);
		Nodes[second_read].add_arc(oe1, 0);
		Nodes[second_read].add_arc(oe2, 0);
    }
    else if( (second_end - second_start) >= Nodes[second_read].len )		// second read is contained
    {
    	assert(second_start == 0);
    	assert(second_end == Nodes[second_read].len);
		Nodes[second_read].flag |= CONTAINED;

		if(!keep)
			return;

        if(complement)
        {
        	oe1 = new OvlEdge(first_read, second_read, first_end - Nodes[first_read].len, first_start, OUTIE|NON_PROPER, OUTIE);
			oe2 = new OvlEdge(first_read, second_read, 0 - first_start, Nodes[first_read].len - first_end, INNIE|NON_PROPER, INNIE);
        }
        else
        {
        	oe1 = new OvlEdge(first_read, second_read, first_end - Nodes[first_read].len, first_start, OUTIE|NON_PROPER, INNIE);
        	oe2 = new OvlEdge(first_read, second_read, 0 - first_start, Nodes[first_read].len - first_end, INNIE|NON_PROPER, OUTIE);
        }

		Nodes[first_read].add_arc(oe1, 1);
		Nodes[first_read].add_arc(oe2, 1);
		Nodes[second_read].add_arc(oe1, 0);
		Nodes[second_read].add_arc(oe2, 0);
    }
    else	// proper overlap
    {
    	if(first_read == second_read)	// does not allow a self-edge!
			return;

		if(second_start == 0)	// if the arc is from the first read to the second
		{
			if(complement)
				oe1 = new OvlEdge(first_read, second_read, Nodes[second_read].len - second_end, first_start, OUTIE, OUTIE);
			else
				oe1 = new OvlEdge(first_read, second_read, Nodes[second_read].len - second_end, first_start, OUTIE, INNIE);

			Nodes[first_read].add_arc(oe1, 1);
			Nodes[second_read].add_arc(oe1, 0);

			assert(Nodes[second_read].len - second_end > 0);
			assert(first_start > 0);

			assert(Nodes[second_read].len - second_end < Nodes[second_read].len);
			assert(first_start < Nodes[first_read].len);

		}
		else if(first_start == 0)	// the arc is from second to first
		{
			if(complement)
				oe1 = new OvlEdge(first_read, second_read, second_start, Nodes[first_read].len - first_end, INNIE, INNIE);
			else
				oe1 = new OvlEdge(first_read, second_read, second_start, Nodes[first_read].len - first_end, INNIE, OUTIE);

			Nodes[first_read].add_arc(oe1, 1);
			Nodes[second_read].add_arc(oe1, 0);

			assert(Nodes[first_read].len - first_end > 0);
			assert(second_start > 0);

			assert( Nodes[first_read].len - first_end < Nodes[first_read].len);
			assert(second_start < Nodes[second_read].len);
		}
		else
		{
			throw "OvlGraph: Overlaps are corrupted! Please check the file!";
		}
    }
}

int OvlGraph::venom()
{
	int * neighbourCount = new int[numNodes+1];
	int * neighbourCache = new int[numNodes+1];

	for(int A=1; A<=numNodes; A++)
		neighbourCache[A] = 0;

	for(int A=1; A<=numNodes; A++)
	{
		for(int d=0; d<2; d++)
		{
			Direction dir1 = (d==0) ? INNIE : OUTIE;

			LinkedIter<OvlEdge> iter_A(Nodes[A].get_edge_list(dir1));

			for(OvlEdge * oe_A = iter_A.get_first(); oe_A != NULL; oe_A = iter_A.get_next())
			{
				int i = iter_A.get_flg();
				int B = oe_A->toRead[i];
				int Astart = oe_A->len[i];
				int Bstart = oe_A->len[(i+1)%2];

				int countA = 0;
				LinkedIter<OvlEdge> iter_A2(Nodes[A].get_edge_list(dir1));
				for(OvlEdge * oe_A2 = iter_A2.get_first(); oe_A2 != NULL; oe_A2 = iter_A2.get_next())
				{
					int j = iter_A2.get_flg();
					int C = oe_A2->toRead[j];

					if(!oe_A2->flagged(REMOVABLE) && oe_A2->len[j] < Astart)
					{
						neighbourCount[countA++] = C;
						neighbourCache[C] = 1;
					}
				}

				Direction dir2 = (oe_A->arrow[i] & OUTIE) ? INNIE : OUTIE;

				bool toremove = true;

				int countB = 0;
				LinkedIter<OvlEdge> iter_B(Nodes[B].get_edge_list(dir2));
				for(OvlEdge * oe_B = iter_B.get_first(); oe_B != NULL; oe_B = iter_B.get_next())
				{
					int j = iter_B.get_flg();
					int C = oe_B->toRead[j];

					if(!oe_B->flagged(REMOVABLE) && oe_B->len[j] < Bstart)
					{
						countB++;
						if(neighbourCache[C] == 1)	// we have a hit
							toremove = false;
					}
				}

				if(countA>0 && countB>0 && toremove)
					oe_A->flag(REMOVABLE);

				for(int k=0; k<countA; k++)
					neighbourCache[ neighbourCount[k] ] = 0;
			}
		}
	}

	delete [] neighbourCount;
	delete [] neighbourCache;

	int counter = 0;
	// Nothing is actually removed till now. Next we will remove everything marked to be removed
	for(int A=1; A<=numNodes; A++)
		counter += Nodes[A].remove_marked(REMOVABLE);

	return counter;
}

// remove mate edges that are shorter than the given cutoff
void OvlGraph::remove_short_mate_overlaps(int cutoff)
{
	for(int x=1; x<=numPairs; x++)
	{
		int B = 2*x-1;
		LinkedIter<MateEdge> iter( Pairs[x].outArcs );

		for(int d=0; d<2; d++)
		{
			if(d)
			{
				B = 2*x;
				iter.change_list(Pairs[x].inArcs);
			}

			int numLongs = 0;
			for(MateEdge * me = iter.get_first(); me != NULL; me = iter.get_next())
			{
				int i = iter.get_flg();
				int y = me->toRead[i];
				int C = (me->arrow[i] & INNIE) ? 2*y : 2*y-1;

				if(Nodes[B].len + me->len[2] + Nodes[C].len >= cutoff)
					numLongs++;
			}

			if(numLongs > 0)
			{
				for(MateEdge * me = iter.get_first(); me != NULL; me = iter.get_next())
				{
					int i = iter.get_flg();
					int y = me->toRead[i];
					int C = (me->arrow[i] & INNIE) ? 2*y : 2*y-1;

					if(Nodes[B].len + me->len[2] + Nodes[C].len < cutoff)
						me->flag(REMOVABLE, i);
				}
			}
		}
	}
	prune_graph<MatePair>(Pairs, numPairs, REMOVABLE);
}

void OvlGraph::remove_short_read_overlaps(int cutoff)
{
	cutoff *= 2;
	for(int x=1; x<=numNodes; x++)
	{
		LinkedIter<OvlEdge> iter( Nodes[x].outArcs );

		for(int d=0; d<2; d++)
		{
			if(d) iter.change_list(Nodes[x].inArcs);

			int numLongs = 0;
			for(OvlEdge * me = iter.get_first(); me != NULL; me = iter.get_next())
			{
				int i = iter.get_flg();
				int y = me->toRead[i];

				if(Nodes[x].len - me->len[0] + Nodes[y].len - me->len[1] >= cutoff)
					numLongs++;
			}

			if(numLongs> 0)
			{
				for(OvlEdge * me = iter.get_first(); me != NULL; me = iter.get_next())
				{
					int i = iter.get_flg();
					int y = me->toRead[i];
					if(Nodes[x].len - me->len[0] + Nodes[y].len - me->len[1] < cutoff)
						me->flag(REMOVABLE, i);
				}
			}
		}
	}
	prune_graph<ReadNode>(Nodes, numNodes, REMOVABLE);	// it will remove an edge only if both arrows are flagged
}

// mark nodes that seem to belong to unique regions
void OvlGraph::mark_single_nodes(double arrival, double cutoff)
{
	double ln2 = 0.69;
	for(int x=1; x<=numNodes; x++)
	{
		LinkedIter<PathEdge> iter(Reads[x].outArcs);
		for(int d=0; d<2; d++)
		{
			if(d)	iter.change_list(Reads[x].inArcs);
			for(PathEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
			{
				if((e_ptr->len/arrival) - e_ptr->path_size()*ln2 > cutoff)
				{
					LinkedIter<EdgeTuple> iter((e_ptr->path[0]));
					for(EdgeTuple * et = iter.get_first(); et != NULL; et = iter.get_next())
						Nodes[et->read].flag |= SINGLE;
				}
			}
		}
	}
}

// adjust the insert size statistics
void OvlGraph::calibrate_libraries(char * filename)
{
	long int * libMean = new long int[numLibs+1];
	int * libSupport = new int[numLibs+1];
	int * libCov = new int[numLibs+1];

	int histSize = 100000;
	int ** histogram = new int*[numLibs+1];

	for(int k=0; k<=numLibs; k++)
	{
		libMean[k] = 0;
		libSupport[k] = 0;
		libCov[k] = 0;

		histogram[k] = new int[histSize+1];
		for(int d=0; d<=histSize; d++)
			histogram[k][d] = 0;
	}

	int * offsets = new int[numNodes+1];
	bool * orients = new bool[numNodes+1];

	int end = 2*numPairs;

	for(int x=0; x<=end; x++)
	{
		offsets[x] = -1;
		orients[x] = true;
		if(x%2==1)
		{
			int t = Pairs[(x+1)/2].lib;
			if(Nodes[x].get_degree() > 0 || Nodes[x+1].get_degree() > 0)
				libCov[t] += 1;
		}
	}

	for(int x=1; x<=end; x++)
	{
		LinkedIter<PathEdge> iter(Reads[x].outArcs);
		for(int d=0; d<2; d++)
		{
			if(d)	iter.change_list(Reads[x].inArcs);
			for(PathEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
			{
				int pos = 0;
				LinkedIter<EdgeTuple> iter(e_ptr->path[0]);
				for(EdgeTuple * et = iter.get_first(); et != NULL; et = iter.get_next())
				{
					int A = et->read;
					int B = (A%2==0) ? A - 1 : A + 1;

					pos += et->offset;
					offsets[A] = pos;
					orients[A] = (et->orientation & INNIE);

					if(!orients[A] && offsets[B] >= 0 && orients[B])
					{
						int t = Pairs[ (A+B+1)/4 ].lib;
						int distance = pos + Nodes[A].len - offsets[B];
						if(distance < histSize)
						{
							histogram[t][distance] += 1;
							libSupport[t] += 1;
							libMean[t] += distance;
						}
						else
							histogram[t][histSize] += 1;	// it is clearly an outlier
					}
				}

				for(EdgeTuple * et = iter.get_first(); et != NULL; et = iter.get_next())	// clean up
					offsets[et->read] = -1;
			}
		}
	}

	delete [] offsets;
	delete [] orients;

	for(int k=1; k<=numLibs; k++)
	{
		double confidence = (double)libSupport[k]/libCov[k];

		cout << "--- Statistics for library #" << k;
		cout << " (support: " << libSupport[k] <<", confidence: "<< confidence << ") --- " << endl;
		if(confidence < 0.5)	// confidence is too low
		{
			cout << "WARNING: confidence is too low, library statistics will not be altered!" << endl;
			continue;
		}

		libMean[k] = int(libMean[k]/(double)libSupport[k]);
		bootstrap_library(libMean, libSupport, histogram, histSize, Libs, k, true);
	}
	
	string s(filename);		
	s += ".hist";
	char * fname = new char[s.size()+1];
	strcpy(fname, s.c_str());
	ofstream fh;
	open_n_check(fh, fname);
	
	for(int d=0; d<histSize; d++)
	{
		for(int k=1; k<=numLibs; k++)
			fh << histogram[k][d] << " ";
		fh << "\n";
	}
	check_n_close(fh);
	delete [] fname;

	for(int k=1; k<=numLibs; k++)
		delete [] histogram[k];

	delete [] histogram;
	delete [] libMean;
	delete [] libSupport;
	delete [] libCov;
}

// create an artificial bridge between matepairs that have ends on the same path
void OvlGraph::bridge_matepairs()
{
	for(int x=1; x<=numPairs; x++)
	{
		if(Pairs[x].get_degree() == 1)
		{
			int B = (Pairs[x].get_out_degree() == 1) ? 2*x : 2*x-1;
			Direction dir = OUTIE;
			Direction dir2 = (B%2==0) ? INNIE : OUTIE;

			int dist = 0;

			if(Nodes[B].isChild())
			{
				if(!Nodes[B].same) dir = INNIE;
				dist = Nodes[B].offset;
				B = Nodes[B].parent;
			}

			if(Nodes[B].is_middle())
			{
				while(true)
				{
					LinkedIter<OvlEdge> iter( Nodes[B].get_edge_list(dir) );
					OvlEdge * e_ptr = iter.get_first();

					if(e_ptr == NULL)
						break;

					int k = iter.get_flg();
					int C = e_ptr->toRead[k];

					dist += e_ptr->len[k];

					int y = (C%2==0) ? C/2 : (C+1)/2;

					if(x != y && Nodes[C].is_middle() && B != C)
					{
						dir = (e_ptr->arrow[k] & INNIE) ? INNIE : OUTIE;
						B = C;

						if(y <= numPairs && dir == INNIE && Pairs[y].get_degree() > 0 && !Pairs[y].is_middle())
						{
							Direction dir3 = (C%2==0) ? INNIE : OUTIE;

							dir2 |= FAKE;
							dir3 |= FAKE;

							MateEdge * me = new MateEdge(x, y, dir2, dir3, dist + Nodes[C].len, dist + Nodes[C].len, dist + Nodes[C].len);

							Pairs[x].add_arc(me, 1);
							Pairs[y].add_arc(me, 0);

							break;
						}
					}
					else
						break;
				}
			}
		}
	}
}

int OvlGraph::print_mate_contigs(char * filename)
{
	ofstream fh;
	open_n_check(fh, filename);

	int problem = 0, gaps = 0;

	QStruct * que = new QStruct[numNodes+1];
	PStruct * ptrs = new PStruct[numNodes+1];
	SDNode * seen = new SDNode[numNodes+1];
	int * in_parents = new int[numNodes+1];
	int * out_parents = new int[numNodes+1];
	int * reported = new int[numNodes+1];

	for(int x=0; x<=numNodes; x++)
	{
		ptrs[x].x[0] = -1;
		ptrs[x].x[1] = -1;
		ptrs[x].x[2] = -1;

		que[x].readID = 0;
		que[x].dir = -1;
		que[x].val = -1;

		seen[x].target = 0;
		seen[x].in = -1;
		seen[x].out = -1;

		reported[x] = 0;
	}

	bool * labels = new bool[numPairs+1];
	for(int x=1; x<=numPairs; x++)
		labels[x] = true;

	for(int x=1; x<=numPairs; x++)
	{
		int A = 2*x;
		LinkedIter<PathEdge> iter(Mates[x].outArcs);

		for(int d=0; d<2; d++)
		{
			if(d)
			{
				A = 2*x-1;
				iter.change_list(Mates[x].inArcs);
			}

			for(PathEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
			{
				int i = iter.get_flg();
				int y = e_ptr->toRead[i];

				if( e_ptr->path[i].get_size() < 2 || (x > y || (x==y && i==1)) )
					continue;

				int next, before = A, prev = A;
				LinkedIter<EdgeTuple> et_iter(e_ptr->path[i]);

				bool first = true;
				bool isFake = true;
				bool isRepeat = true;

				for(EdgeTuple * et = et_iter.get_first(); et != NULL; et = et_iter.get_next())
				{
					next = (et->orientation & INNIE) ? 2*et->read : 2*et->read - 1;

					if(!first || labels[x])
					{
						labels[x] = false;
						if(et->orientation & FAKE)
						{
							int prev2 = (prev%2==0) ? prev-1 : prev+1;
							int z = (prev%2==0) ? prev/2 : (prev+1)/2;

							int L = Pairs[z].lib;

							gaps += close_gap(fh, prev, prev2, Libs[L].maxDist,
								que, seen, ptrs, in_parents, out_parents, ++problem, INNIE, OUTIE, reported, isRepeat);
							gaps += close_gap(fh, prev2, next, int((et->offset + Nodes[prev2].len + Nodes[next].len)),
								que, seen, ptrs, in_parents, out_parents, ++problem, OUTIE, INNIE, reported, isRepeat);

							isFake = true;
						}
						else
						{
							gaps += close_gap(fh, prev, next, int((et->offset + Nodes[prev].len + Nodes[next].len)),
								que, seen, ptrs, in_parents, out_parents, ++problem, INNIE, INNIE, reported, isRepeat);

							isFake = false;
						}
					}

					before = prev;
					prev = next;
					first = false;
				}

				int z = (prev%2==0) ? prev/2 : (prev+1)/2;
				assert(z == y);

				if(labels[z])
				{
					labels[z] = false;
					int L = Pairs[z].lib;

					int prev2 = (prev%2==0) ? prev-1 : prev+1;
					int after = (before%2==0) ? before-1 : before+1;

					if(isFake)
					{
						gaps += close_gap(fh, prev, prev2, Libs[L].maxDist,
							que, seen, ptrs, in_parents, out_parents, ++problem, INNIE, OUTIE, reported, isRepeat);
					}
					else
					{
						gaps += close_gap(fh, prev, after, Libs[L].maxDist,
							que, seen, ptrs, in_parents, out_parents, ++problem, INNIE, OUTIE, reported, isRepeat);

						gaps += close_gap(fh, after, prev2, Libs[L].maxDist,
							que, seen, ptrs, in_parents, out_parents, ++problem, OUTIE, OUTIE, reported, isRepeat);
					}

					int B = prev2;
					Direction dir = OUTIE;

					if(Nodes[B].isChild())
					{
						if(!Nodes[B].same) dir = INNIE;
						B = Nodes[B].parent;
					}

					if(!isRepeat && Nodes[B].is_middle())	// if it's on a long chain, then report until the end of the chain
					{
						while(true)
						{
							LinkedIter<OvlEdge> iter( Nodes[B].get_edge_list(dir) );
							OvlEdge * e_ptr = iter.get_first();

							if(e_ptr == NULL)
								break;

							int k = iter.get_flg();
							int C = e_ptr->toRead[k];

							if(Nodes[C].is_middle() && !reported[C] && B != C)
							{
								dir = (e_ptr->arrow[k] & INNIE) ? INNIE : OUTIE;
								fh << e_ptr->len[k] << " " << dir << " " << C << "\n";
								reported[C] += 1;
								B = C;
							}
							else
								break;
						}
					}
				}

				if(isRepeat)
					fh << "-1 -1" << "\n";
				else
					fh << "-1 0" << "\n";

			}
		}
	}

	delete [] ptrs;
	delete [] que;
	delete [] seen;
	delete [] in_parents;
	delete [] out_parents;
	delete [] labels;

	for(int i=1; i<=numNodes; i++)
		Reads[i].print_edges(fh, reported, 5);

	check_n_close(fh);

	delete [] reported;
	return gaps;
}

int OvlGraph::check_consistency()
{
	int problem = 0, counter = 0;

	QStruct * que = new QStruct[numNodes+1];
	PStruct * ptrs = new PStruct[numNodes+1];
	SDNode * seen = new SDNode[numNodes+1];

	int * roots = new int[numNodes+1];

	for(int x=0; x<=numNodes; x++)
	{
		ptrs[x].x[0] = -1;
		ptrs[x].x[1] = -1;
		ptrs[x].x[2] = -1;

		que[x].readID = 0;
		que[x].dir = -1;
		que[x].val = -1;

		seen[x].target = 0;
		seen[x].in = -1;
		seen[x].out = -1;

		roots[x] = x;
	}

	int bads = 0;

	for(int x=1; x<=numPairs; x++)
	{
		if(Pairs[x].get_in_degree() == 0 || Pairs[x].get_out_degree() == 0)
			continue;

		for(int d=0; d<2; d++)
		{
			LinkedIter<MateEdge> iter(Pairs[x].inArcs);
			if(d) iter.change_list(Pairs[x].outArcs);

			for(MateEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
			{
				int i = iter.get_flg();
				int y = e_ptr->toRead[i];

				bool remove = true;

				LinkedIter<MateEdge> iter2(Pairs[x].outArcs);
				if(d) iter2.change_list(Pairs[x].inArcs);

				for(MateEdge * e_ptr2 = iter2.get_first(); e_ptr2 != NULL; e_ptr2 = iter2.get_next())
				{
					int j = iter2.get_flg();
					int z = e_ptr2->toRead[j];

					if( !Pairs[y].isNeighbour(z, e_ptr->arrow[i] & (INNIE|OUTIE), e_ptr2->arrow[j] & (OUTIE|INNIE) ) )
					{
						int source = (e_ptr->arrow[i] & OUTIE) ? 2*y-1 : 2*y;
						int target = (e_ptr2->arrow[j] & OUTIE) ? 2*z-1 : 2*z;
						int maxdist = Pairs[x].len - e_ptr->len[2] - e_ptr2->len[2];

						if( isConsistent(source, target, maxdist, que, seen, ptrs, ++problem, OUTIE, INNIE) ||
							isConsistent(source, target, Pairs[x].len, que, seen, ptrs, ++problem, INNIE, OUTIE) )
						{
							remove = false;

							int rz = z;
							while(roots[rz] != rz) rz = roots[rz];

							int ry = y;
							while(roots[ry] != ry) ry = roots[ry];

							roots[rz] = ry;		// does not matter which direction
						}
						else
							bads++;
					}
					else
					{
						remove = false;

						int rz = z;
						while(roots[rz] != rz) rz = roots[rz];

						int ry = y;
						while(roots[ry] != ry) ry = roots[ry];

						roots[rz] = ry;		// does not matter which direction
					}
				}

				if(iter2.get_size() > 0 && remove)
				{
					e_ptr->flag(REMOVABLE);
					Pairs[y].remove_marked(REMOVABLE);
					roots[y] = y;
					counter++;
				}
			}
			Pairs[x].remove_marked(REMOVABLE);
		}

		if(Pairs[x].get_in_degree() > 0 && Pairs[x].get_out_degree() > 0)
		{
			int root = 0;
			int pro = 1;
			int con = 0;

			for(int d=0; d<2; d++)
			{
				LinkedIter<MateEdge> iter(Pairs[x].inArcs);
				if(d) iter.change_list(Pairs[x].outArcs);

				for(MateEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
				{
					int i = iter.get_flg();
					int y = e_ptr->toRead[i];

					int k = y;
					while(roots[k] != k) k = roots[k];

					if(root == 0)
						root = k;
					else if(root != k)
						con++;
					else
						pro++;
				}
			}

			if(con > 1 && pro > 1 && pro != con)
			{
				LinkedIter<MateEdge> iter(Pairs[x].inArcs);
				for(int d=0; d<2; d++)
				{
					if(d) iter.change_list(Pairs[x].outArcs);
					for(MateEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
					{
						int i = iter.get_flg();
						int y = e_ptr->toRead[i];

						int k = y;
						while(roots[k] != k) k = roots[k];

						if( (pro > con && k != root) || (pro < con && k == root) )
						{
							e_ptr->flag(REMOVABLE);
							Pairs[y].remove_marked(REMOVABLE);
							roots[y] = y;
							counter++;
						}
					}
				}
				Pairs[x].remove_marked(REMOVABLE);
			}
		}

		for(int d=0; d<2; d++)
		{
			LinkedIter<MateEdge> iter(Pairs[x].inArcs);
			if(d) iter.change_list(Pairs[x].outArcs);

			for(MateEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
			{
				int i = iter.get_flg();
				int k = e_ptr->toRead[i];
				roots[k] = k;
			}
		}
	}

	delete [] ptrs;
	delete [] que;
	delete [] seen;
	delete [] roots;

	return bads;
}

bool OvlGraph::isConsistent(int source, int target, int maxdist, QStruct * que, SDNode * seen, PStruct * ptrs,
	int problem, Direction sourcedir, Direction targetdir)
{
	int que_size = 1;
	maxdist = int(maxdist * MatePair::FLEX);

	que[0].readID = source;
	que[0].dir = sourcedir;
	que[0].val = 0;

	if(sourcedir == OUTIE)
		source = 0 - source;

	int not_in_que = -1;

	while(que_size > 0)
	{
		int A = que[0].readID;
		Direction dir = que[0].dir;
        int val = que[0].val;

		if(seen[A].target != problem)	// if we encounter A the first time for this instance
		{
			seen[A].target = problem;
			seen[A].in = -1;
			seen[A].out = -1;
		}

		if(dir == OUTIE)
			seen[A].out = val;
		else
			seen[A].in = val;

		if(A == target && dir == targetdir)	// then we're done
		{
			for(int j=0; j<que_size; j++)
				ptrs[ que[j].readID ].x[ que[j].dir ] = not_in_que;	// pop out the rest of the queue

			return true;
		}

		swap_que_elts(0, --que_size, que, ptrs);

		// remove the min
		int j = 0;

		int leftchild = 2*(j+1)-1;
		int rightchild = leftchild+1;

		while( leftchild < que_size )
		{
			int parent_val = que[j].val;
			int leftchild_val = que[leftchild].val;
			int rightchild_val = (rightchild < que_size) ? que[rightchild].val : leftchild_val + 1;

			if(parent_val <= leftchild_val && parent_val <= rightchild_val)
				break;

			if(rightchild_val < leftchild_val)
			{
				swap_que_elts(j, rightchild, que, ptrs);
				j = rightchild;
			}
			else
			{
				swap_que_elts(j, leftchild, que, ptrs);
				j = leftchild;
			}
			leftchild = 2*(j+1)-1;
			rightchild = leftchild + 1;
		}
		ptrs[ que[j].readID ].x[ que[j].dir ] = j;

		ptrs[A].x[dir] = not_in_que;	// pop A out of the queue

		LinkedIter<OvlEdge> iter( Nodes[A].get_edge_list(dir) );

		for(OvlEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
		{
			int k = iter.get_flg();
			int C = e_ptr->toRead[k];

			int stamp = val + e_ptr->len[k];
			Direction new_dir = (e_ptr->arrow[k] & INNIE) ? INNIE : OUTIE;

			int i = -1;

			// if not in que and not been processed for this direction
			if(ptrs[C].x[new_dir] == not_in_que && (seen[C].target != problem || ((new_dir == OUTIE && seen[C].out == -1) || (new_dir == INNIE && seen[C].in == -1))))
			{
				if( (stamp + Nodes[C].len) <= maxdist )
				{
					i = que_size;

					que[i].readID = C;
					que[i].dir = new_dir;
					que[i].val = stamp;

					ptrs[C].x[new_dir] = i;

					que_size++;
				}
			}
			else if(ptrs[C].x[new_dir] != not_in_que)
			{
				i = ptrs[C].x[new_dir];		// get its location in the queue
				if( que[i].val > stamp)		// a better route
					que[i].val = stamp;
				else
					i = -1;
			}

			if(i > -1)
			{
				while(i > 0)
				{
					int p = int((i-1)/2);
					int parent_val = que[p].val;

					if(parent_val <= stamp)
					{	break;	}

					swap_que_elts(i, p, que, ptrs);
					i = p;
				}
				que[i].readID = C;
				que[i].dir = new_dir;
				que[i].val = stamp;
				ptrs[C].x[new_dir] = i;
			}
		}
	}

	return false;
}


// close the gap between the mate nodes
int OvlGraph::close_gap(ofstream & fh, int source, int target, int maxdist, QStruct * que, SDNode * seen, PStruct * ptrs,
	int * in_parents, int * out_parents, int problem, Direction sourcedir, Direction targetdir, int * reported, bool & isRepeat)
{
	int que_size = 1;
	maxdist = int(maxdist * MatePair::FLEX);

	if(source == Nodes[target].parent)	// then just pretend that we processed, nothing to do
		return 0;

	if(Nodes[source].isChild())	// if the source is a contained read
	{
		if(sourcedir == INNIE)
			maxdist += Nodes[source].offset;
		else
			maxdist += (Nodes[Nodes[source].parent].len - Nodes[source].len - Nodes[source].offset);

		if(!Nodes[source].same)
			sourcedir = (sourcedir & INNIE) ? OUTIE : INNIE;
		source = Nodes[source].parent;
	}

	if(Nodes[target].isChild())	// if the target is a contained read
	{
		if(targetdir == INNIE)
			maxdist += (Nodes[Nodes[target].parent].len - Nodes[target].len - Nodes[target].offset);
		else
			maxdist += Nodes[target].offset;

		if(!Nodes[target].same)
			targetdir = (targetdir & INNIE) ? OUTIE : INNIE;
		target = Nodes[target].parent;
	}

	que[0].readID = source;
	que[0].dir = sourcedir;
	que[0].val = 0;

	//reported[source] += 1;
	//reported[target] += 1;

	if(sourcedir == OUTIE)
		source = 0 - source;

	int not_in_que = -1;

	while(que_size > 0)
	{
		int A = que[0].readID;
		Direction dir = que[0].dir;
        int val = que[0].val;

		if(seen[A].target != problem)	// if we encounter A the first time for this instance
		{
			seen[A].target = problem;
			seen[A].in = -1;
			seen[A].out = -1;
		}

		if(dir == OUTIE)
			seen[A].out = val;
		else
			seen[A].in = val;

		if(A == target && dir == targetdir)	// then we're done
		{
			int x = (dir == INNIE) ? A : -A;
			int y;
			int counter = -1;

			for(int j=0; j<que_size; j++)
				ptrs[ que[j].readID ].x[ que[j].dir ] = not_in_que;	// pop out the rest of the queue

			while(x != source)
			{
				que[++counter].readID = abs(x);
				if(x<0)
				{
					y = out_parents[abs(x)];
					que[counter].dir = OUTIE;
					que[counter].val = seen[abs(x)].out;
				}
				else
				{
					y = in_parents[x];
					que[counter].dir = INNIE;
					que[counter].val = seen[x].in;
				}
				que[counter].val = (y>0) ? que[counter].val - seen[y].in : que[counter].val - seen[abs(y)].out ;
				x = y;
			}

			for(int i=counter; i>=0; i--)
			{
				if(que[i].val < 0) que[i].val = 1;
				fh << que[i].val << " " << que[i].dir << " " << que[i].readID << "\n";
				Nodes[que[i].readID].print_contains(fh);
				if(reported[que[i].readID] == 0)
					isRepeat = false;
				reported[que[i].readID] += 1;
			}

			return 0;
		}

		swap_que_elts(0, --que_size, que, ptrs);

		// remove the min
		int j = 0;

		int leftchild = 2*(j+1)-1;
		int rightchild = leftchild+1;

		while( leftchild < que_size )
		{
			int parent_val = que[j].val;
			int leftchild_val = que[leftchild].val;
			int rightchild_val = (rightchild < que_size) ? que[rightchild].val : leftchild_val + 1;

			if(parent_val <= leftchild_val && parent_val <= rightchild_val)
				break;

			if(rightchild_val < leftchild_val)
			{
				swap_que_elts(j, rightchild, que, ptrs);
				j = rightchild;
			}
			else
			{
				swap_que_elts(j, leftchild, que, ptrs);
				j = leftchild;
			}
			leftchild = 2*(j+1)-1;
			rightchild = leftchild + 1;
		}
		ptrs[ que[j].readID ].x[ que[j].dir ] = j;

		ptrs[A].x[dir] = not_in_que;	// pop A out of the queue

		LinkedIter<OvlEdge> iter( Nodes[A].get_edge_list(dir) );

		for(OvlEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
		{
			int k = iter.get_flg();
			int C = e_ptr->toRead[k];

			int stamp = val + e_ptr->len[k];
			Direction new_dir = (e_ptr->arrow[k] & INNIE) ? INNIE : OUTIE;

			int i = -1;

			// if not in que and not been processed for this direction
			if(ptrs[C].x[new_dir] == not_in_que && (seen[C].target != problem || ((new_dir == OUTIE && seen[C].out == -1) || (new_dir == INNIE && seen[C].in == -1))))
			{
				if( (stamp + Nodes[C].len) <= maxdist )
				{
					i = que_size;

					que[i].readID = C;
					que[i].dir = new_dir;
					que[i].val = stamp;
					if(new_dir == INNIE)
						in_parents[C] = (dir == INNIE) ? A : -A;
					else
						out_parents[C] = (dir == INNIE) ? A : -A;

					ptrs[C].x[new_dir] = i;

					que_size++;
				}
			}
			else if(ptrs[C].x[new_dir] != not_in_que)
			{
				i = ptrs[C].x[new_dir];		// get its location in the queue
				if( que[i].val > stamp)		// a better route
				{
					que[i].val = stamp;
					if(new_dir == INNIE)
						in_parents[C] = (dir == INNIE) ? A : -A;
					else
						out_parents[C] = (dir == INNIE) ? A : -A;
				}
				else
					i = -1;
			}

			if(i > -1)
			{
				while(i > 0)
				{
					int p = int((i-1)/2);
					int parent_val = que[p].val;

					if(parent_val <= stamp)
					{	break;	}

					swap_que_elts(i, p, que, ptrs);
					i = p;
				}
				que[i].readID = C;
				que[i].dir = new_dir;
				que[i].val = stamp;
				ptrs[C].x[new_dir] = i;
			}
		}
	}

	fh << "-1 0" << "\n";
	return 1;		// shouldn't reach here
}

// use Dijkstra's algorithm to find the set of nodes between a mate pair
int OvlGraph::verify_distance(int source, int target, int maxdist, QStruct *que, SDNode *seen, PStruct * ptrs,
	                           Direction start_dir, int *seenlist, int &counter, SDNode *seen_before, bool ignore)
{
	int shortest_dist = maxdist;
	int que_size = 1;
	int not_in_que = -1;
	counter = 0;

	que[0].readID = source;		// start with the source
	que[0].dir = start_dir;
	que[0].val = 0;

	if(Nodes[source].isChild())
	{
		que[0].readID = Nodes[source].parent;
		if(!Nodes[source].same)
			que[0].dir = (start_dir & INNIE) ? OUTIE : INNIE;
	}

	int deputy = target;
	if(Nodes[target].isChild())
		deputy = Nodes[target].parent;

	while(que_size > 0)
	{
		int A = que[0].readID;
		Direction dir = que[0].dir;
        EdgeLength val = que[0].val;

		Nodes[A].process(A, seen, seenlist, target, counter, dir, val);

		swap_que_elts(0, --que_size, que, ptrs);

		// remove the min
		int j = 0;

		int leftchild = 2*(j+1)-1;
		int rightchild = leftchild+1;

		while( leftchild < que_size )
		{
			int parent_val = que[j].val;
			int leftchild_val = que[leftchild].val;
			int rightchild_val = (rightchild < que_size) ? que[rightchild].val : leftchild_val + 2;

			if(parent_val <= leftchild_val && parent_val <= rightchild_val)
				break;

			if(rightchild_val < leftchild_val)
			{
				swap_que_elts(j, rightchild, que, ptrs);
				j = rightchild;
			}
			else
			{
				swap_que_elts(j, leftchild, que, ptrs);
				j = leftchild;
			}
			leftchild = 2*(j+1)-1;
			rightchild = leftchild + 1;
		}
		ptrs[ que[j].readID ].x[ que[j].dir ] = j;

		ptrs[A].x[dir] = not_in_que;	// pop A out of the queue

		LinkedIter<OvlEdge> iter(Nodes[A].get_edge_list(dir));

		for(OvlEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
		{
			int k = iter.get_flg();
			int C = e_ptr->toRead[k];
			Direction new_dir = (e_ptr->arrow[k] & INNIE) ? INNIE : OUTIE;

			if(!ignore && seen_before[C].target != source)
				continue;

			int cutoff = (val + e_ptr->len[k] + Nodes[C].len);

			if(C == deputy && cutoff < maxdist)
			{
				if( (Nodes[target].same && e_ptr->arrow[k] & OUTIE) || (!Nodes[target].same && e_ptr->arrow[k] & INNIE) )
				{
					if(cutoff < shortest_dist)
						shortest_dist = cutoff;
				}
			}

			int stamp = val + e_ptr->len[k];
			int i = -1;

			// if not in que and not been processed for this direction
			if(ptrs[C].x[new_dir] == not_in_que && (seen[C].target != target ||
				((new_dir == OUTIE && seen[C].out == -1) || (new_dir == INNIE && seen[C].in == -1))))
			{
				if( (stamp + Nodes[C].len) < maxdist )
				{
					i = que_size;

					que[i].readID = C;
					que[i].dir = new_dir;
					que[i].val = stamp;
					ptrs[C].x[new_dir] = i;

					que_size++;
				}
			}
			else if(ptrs[C].x[new_dir] != not_in_que)
			{
				i = ptrs[C].x[new_dir];		// get its location in the queue
				if( que[i].val > stamp)		// a better route
					que[i].val = stamp;
				else
					i = -1;
			}

			if(i > -1)
			{
				while(i > 0)
				{
					int p = int((i-1)/2);
					int parent_val = que[p].val;

					if(parent_val <= stamp)
					{	break;	}

					swap_que_elts(i, p, que, ptrs);
					i = p;
				}
				que[i].readID = C;
				que[i].dir = new_dir;
				que[i].val = stamp;
				ptrs[C].x[new_dir] = i;
			}

		}
	}
	return shortest_dist;
}

// check if the mate pair starting from A has exactly one path
bool OvlGraph::checkPath(int A, int & pathLength, int * vsted)
{
	int next = A;
	short int dir = INNIE;	// assume we are coming from an in-arrow

	while(next > 0)
	{
		LinkedIter<OvlEdge> iter( Nodes[next].get_edge_list(dir) );
		next = -1;

		for(OvlEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
		{
			short int k = iter.get_flg();
			int C = e_ptr->toRead[k];

			if(vsted[C] > 1)	// than we have visited C before
				return false;

			if(vsted[C] == 1)
			{
				if(next > 0)	// then we have a conflict (i.e a fork)
					return false;

				vsted[C] += 1;

				next = C;		// otherwise next in line is C
				dir = (e_ptr->arrow[k] & INNIE) ? INNIE : OUTIE;
			}
		}
	}

	return true;	// if we didn't return false yet then it means the path is single
}

// compute the mate pair graph
void OvlGraph::build_pair_graph()
{
	int st = 1;
	int end = 2*numPairs;

	QStruct * que = new QStruct[numNodes+1];
	PStruct * ptrs = new PStruct[numNodes+1];
	int * list_in = new int[numNodes+1];
	int * list_out = new int[numNodes+1];

	SDNode * seen_A_out = new SDNode[numNodes+1];
	SDNode * seen_A_in = new SDNode[numNodes+1];
	SDNode * seen_B_out = new SDNode[numNodes+1];
	SDNode * seen_B_in = new SDNode[numNodes+1];

	for(int x=0; x<=numNodes; x++)
	{
		ptrs[x].x[0] = -1;
		ptrs[x].x[1] = -1;
		ptrs[x].x[2] = -1;

		que[x].readID = 0;
		que[x].dir = -1;
		que[x].val = -1;

		seen_A_in[x].target = 0;
		seen_A_in[x].in = -1;
		seen_A_in[x].out = -1;

		seen_A_out[x].target = 0;
		seen_A_out[x].in = -1;
		seen_A_out[x].out = -1;

		seen_B_in[x].target = 0;
		seen_B_in[x].in = -1;
		seen_B_in[x].out = -1;

		seen_B_out[x].target = 0;
		seen_B_out[x].in = -1;
		seen_B_out[x].out = -1;
	}

	int found = 0;
	int not_found = 0;
	int dummy;
	int out_counter;
	int in_counter;

	for(int k=1; k<=numLibs; k++)
		Libs[k].maxDist = int((Libs[k].insertSize + Libs[k].sdwidth * Libs[k].deviation) * MatePair::FLEX);

	for(int x=st; x<end; x+=2)
	{
		int A = x;
		int B = x+1;
		int z = int(B/2);

		int pA = Nodes[A].parent;
		int pB = Nodes[B].parent;

		int L = Pairs[z].lib;

		int BA_dist = Libs[L].maxDist;
		Pairs[z].len = BA_dist;

		if( !Pairs[z].isEligible() || (Nodes[A].get_out_degree() == 0 && Nodes[pA].get_degree() == 0) ||
			(Nodes[B].get_out_degree() == 0 && Nodes[pB].get_degree() == 0) )
			continue;

		BA_dist += Nodes[A].offset;	// if they are not contained, the offsets will be zero
		BA_dist += Nodes[B].offset;

		if( (pA%2==0 && pB==(pA-1) && pA/2 <= numPairs && Pairs[pA/2].isEligible()) ||
			(pB%2==0 && pA==(pB-1) && pB/2 <= numPairs && Pairs[pB/2].isEligible()) )
		{
			Pairs[z].flag |= CONTAINED;		// then z is trivially contained
			continue;
		}

		int shortest_dist = verify_distance(A, B, BA_dist, que, seen_A_in, ptrs, INNIE, list_in, dummy, NULL, true);
		if( shortest_dist == BA_dist )		// then no point in trying
		{
			not_found++;
			continue;
		}

		found++;
		Pairs[z].flag |= FOUND;

		verify_distance(A, B, Libs[L].insertSize, que, seen_A_out, ptrs, OUTIE, list_out, dummy, NULL, true);
		verify_distance(B, A, Libs[L].insertSize, que, seen_B_out, ptrs, OUTIE, list_out, out_counter, NULL, true);
		verify_distance(B, A, BA_dist, que, seen_B_in, ptrs, INNIE, list_in, in_counter, seen_A_in, false);

		bool flagged = false;
		for(int k=0; k<out_counter; k++)	// first we're going to check if AB is contained by another mate pair
		{
			int C = list_out[k];
			if(C != A && C != B)
			{
				int D = (C%2==0) ? C-1 : C+1;
				int y = (D%2==0) ? int(D/2) : int(C/2);

				if(y > numPairs || !Pairs[y].isEligible())
					continue;

				int DC_dist = Libs[Pairs[y].lib].maxDist;

				if( (seen_A_out[D].target == B && seen_A_out[D].out > -1 && seen_B_out[C].out > -1) &&
					(seen_A_out[D].out + seen_B_out[C].out + shortest_dist - Nodes[A].len - Nodes[B].len < DC_dist) )
				{
					Pairs[z].flag |= CONTAINED;
					flagged = true;
					break;
				}
			}
		}

		if(flagged) continue;

		for(int k=0; k<in_counter; k++)
		{
			int C = list_in[k];

			int D = (C%2==0) ? C-1 : C+1;
			int y = (D%2==0) ? int(D/2) : int(C/2);

			if( y <= z || y > numPairs || !Pairs[y].isEligible())
				continue;

			int DC_dist = Libs[Pairs[y].lib].maxDist;
			Direction dir = (D%2==0) ? OUTIE : INNIE;

			if(seen_A_in[C].target == B && seen_A_in[C].out > -1 && seen_B_in[C].in > -1)
			{
				int dist = seen_A_in[C].out + Nodes[C].len + seen_B_in[C].in;

				if( dist < BA_dist &&  seen_A_out[D].target == B && seen_A_out[D].out > -1 &&
					Nodes[C].len + seen_A_in[C].out + seen_A_out[D].out - Nodes[A].len + Nodes[D].len < DC_dist )
				{
					int BinCin = seen_B_in[C].in - (Nodes[B].len + Nodes[B].offset);
					int AinCout = seen_A_in[C].out - (Nodes[A].len + Nodes[A].offset);
					int AoutDout = seen_A_out[D].out - (Nodes[pA].len - Nodes[A].offset);

					MateEdge * me = new MateEdge(z, y, OUTIE, dir, BinCin, AinCout, AoutDout);
					Pairs[z].add_arc(me, 1);
					Pairs[y].add_arc(me, 0);
				}
			}

			if(seen_A_in[C].target == B && seen_A_in[C].in > -1 && seen_B_in[C].out > -1)
			{
				int dist = seen_A_in[C].in + Nodes[C].len + seen_B_in[C].out;

				if( dist < BA_dist && seen_B_out[D].target == A && seen_B_out[D].out > -1 &&
					Nodes[C].len + seen_B_in[C].out + seen_B_out[D].out - Nodes[B].len + Nodes[D].len < DC_dist )
				{
					int AinCin = seen_A_in[C].in - (Nodes[A].len + Nodes[A].offset);
					int BinCout = seen_B_in[C].out - (Nodes[B].len + Nodes[B].offset);
					int BoutDout = seen_B_out[D].out - (Nodes[pB].len - Nodes[B].offset);

					MateEdge * me = new MateEdge(z, y, INNIE, dir, AinCin, BinCout, BoutDout);
					Pairs[z].add_arc(me, 1);
					Pairs[y].add_arc(me, 0);
				}
			}
		}
	}

	delete [] ptrs;
	delete [] que;
	delete [] seen_A_in;
	delete [] seen_B_in;
	delete [] seen_A_out;
	delete [] seen_B_out;
	delete [] list_in;
	delete [] list_out;

	cout << "Total number of mate nodes: " << found << endl;
	cout << "Total number of unhappy mates: " << not_found << endl;
}

int OvlGraph::remove_contained_reads_alt()
{
	int num = numNodes;

	OvlEdge ** in_visited = new OvlEdge*[num+1];
	int * in_flags = new int[num+1];

	OvlEdge ** out_visited = new OvlEdge*[num+1];
	int * out_flags = new int[num+1];
	bool * winner = new bool[num+1];

	for(int A=0; A<=num; A++)
	{
		in_visited[A] = NULL;
		out_visited[A] = NULL;
		in_flags[A] = -1;
		out_flags[A] = -1;
		winner[A] = false;
	}

	int removd = 0;

	for(int A=1; A<=num; A++)
	{
		if( !Nodes[A].flag & CONTAINED )	// if node is not contained
			continue;

		bool remove = false;
		LinkedIter<OvlEdge> iter_A(Nodes[A].outArcs);
		for(int d=0; d<2; d++)
		{
			int neighbours = 0;
			if(d)	iter_A.change_list(Nodes[A].inArcs);

			// Step1: Mark the proper neighbours as visited
			OvlEdge * oe_A;
			for(oe_A = iter_A.get_first(); oe_A != NULL; oe_A = iter_A.get_next())
			{
				int i = iter_A.get_flg();
				if( !(oe_A->arrow[i] & NON_PROPER) )
				{
					neighbours++;
					if(oe_A->arrow[i] & INNIE)
					{
						in_visited[ oe_A->toRead[i] ] = oe_A;
						in_flags[ oe_A->toRead[i] ] = i;
					}
					else
					{
						out_visited[ oe_A->toRead[i] ] = oe_A;
						out_flags[ oe_A->toRead[i] ] = i;
					}
				}
			}

			// Step2: Traverse all the parents and see if any contains the child's neighbourhood
			for(oe_A = iter_A.get_first(); oe_A != NULL; oe_A = iter_A.get_next())
			{
				int i = iter_A.get_flg();
				if(oe_A->arrow[i] & NON_PROPER)
				{
					int match = 0;
					int B = oe_A->toRead[i];		// record the information of node B

					LinkedIter<OvlEdge> iter_B(Nodes[B].inArcs);
					if(oe_A->arrow[i] & INNIE)	// if coming inside...
						iter_B.change_list(Nodes[B].outArcs);

					for(OvlEdge * oe_B = iter_B.get_first(); oe_B != NULL; oe_B = iter_B.get_next())
					{
						int j = iter_B.get_flg();
						int C = oe_B->toRead[j];

						if( (in_visited[C] != NULL) && (oe_B->arrow[j] & INNIE) )	// strangely this process is very similar to transitive reduction!
						{
							OvlEdge * oe_C = in_visited[C];
							int k = in_flags[C];

							int fuzz = int((oe_B->len[j] + oe_A->len[i] + oe_C->len[k]) * (MatePair::FLEX-1));
							int diff = oe_B->len[j] + oe_A->len[i] - oe_C->len[k];
							if( abs(diff) <= fuzz )
								match++;
						}
						if( (out_visited[C] != NULL) && (oe_B->arrow[j] & OUTIE) )
						{
							OvlEdge * oe_C = out_visited[C];
							int k = out_flags[C];

							int fuzz = int((oe_B->len[j] + oe_A->len[i] + oe_C->len[k]) * (MatePair::FLEX-1));
							int diff = oe_B->len[j] + oe_A->len[i] - oe_C->len[k];
							if( abs(diff) <= fuzz )
								match++;
						}
					}

					if(match == neighbours)	// bingo!
					{
						if(d && winner[B])	// we can remove the node!
						{
							remove = true;
							Nodes[A].parent = B;
							Nodes[A].same = (oe_A->arrow[i] & INNIE) ? false : true;
							Nodes[A].offset = abs(oe_A->len[rev(i)]);
							Nodes[B].add_containee(A, Nodes[A].same, abs(oe_A->len[rev(i)]), abs(oe_A->len[i]));
							break;
						}
						else
							winner[B] = true;
					}
				}
			}

			// Step 3: Clean up
			for(oe_A = iter_A.get_first(); oe_A != NULL; oe_A = iter_A.get_next())
			{
				int i = iter_A.get_flg();
				int B = oe_A->toRead[i];
				in_visited[ B ] = NULL;
				in_flags[ B ] = -1;
				out_visited[ B ] = NULL;
				out_flags[ B ] = -1;

				if(d)	{ winner[B] = false; }
			}
		}

		remove_node<ReadNode, OvlEdge>(Nodes, A);
		removd++;
	}

	delete [] winner;
	delete [] in_visited;
	delete [] out_visited;
	delete [] in_flags;
	delete [] out_flags;

	return removd;
}
