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

#ifndef OVL_GRAPH_TEMPL_H
#define OVL_GRAPH_TEMPL_H

#include <iostream>
#include <cmath>

#include "PathNode.h"
#include "ReadNode.h"
#include "MatePair.h"
#include "HapUtils.h"
#include "DnaTempl.h"

using namespace std;

// Check if the edge is suitable for reduction (templated to fool the compiler)
template <class Arc>
void check_read_reduction(Arc * oe_B, Arc * oe_A, Arc * oe_C, int i, int j, int k)
{
	int fuzz = int((oe_B->len[j] + oe_A->len[i] + oe_C->len[k]) * (MatePair::FLEX-1));
	int diff = oe_B->len[j] + oe_A->len[i] - oe_C->len[k];
	if( abs(diff) <= fuzz )
	{
		oe_C->flag(REMOVABLE);		// flag as to be removed
		oe_A->flag(SUPPORTED, i);
		oe_B->flag(SUPPORTED, rev(j));
	}
}

// Check if the edge is suitable for reduction (templated to fool the compiler)
template <class Vertex, class Arc>
void check_mate_reduction(Arc * oe_B, Arc * oe_A, Arc * oe_C, int i, int j, int d,
	ReadNode * Nodes, Vertex * Pairs, Direction orient, Direction orient2)
{
	int A = oe_A->toRead[rev(i)];
	int B = oe_B->toRead[rev(j)];
	int C = oe_B->toRead[j];

	bool to_remove = false;

	if(orient & INNIE && orient2 & INNIE)
	{
		if( oe_A->len[i] + Nodes[B*2].len + oe_B->len[j] + Nodes[C*2].len + oe_C->len[2] <= Pairs[A].len)
		{
			if(d)	// processing inArcs
			{
				if( oe_B->len[j] + Nodes[C*2].len + oe_C->len[2] + Nodes[A*2].len + oe_A->len[rev(i)] <= Pairs[B].len )
				{
					if( oe_C->len[2] + Nodes[A*2].len + oe_A->len[rev(i)] + Nodes[B*2-1].len + oe_B->len[rev(j)] <= Pairs[C].len )
						to_remove = true;
				}
			}
			else	// processing outArcs
			{
				if( oe_B->len[j] + Nodes[C*2].len + oe_C->len[2] + Nodes[A*2-1].len + oe_A->len[rev(i)] <= Pairs[B].len )
				{
					if( oe_C->len[2] + Nodes[A*2-1].len + oe_A->len[rev(i)] + Nodes[B*2-1].len + oe_B->len[rev(j)] <= Pairs[C].len )
						to_remove = true;
				}
			}
		}
	}
	else if(orient & INNIE && orient2 & OUTIE)
	{
		if( oe_A->len[i] + Nodes[B*2].len + oe_B->len[j] + Nodes[C*2-1].len + oe_C->len[2] <= Pairs[A].len)
		{
			if(d)	// processing inArcs
			{
				if( oe_B->len[j] + Nodes[C*2-1].len + oe_C->len[2] + Nodes[A*2].len + oe_A->len[rev(i)] <= Pairs[B].len )
				{
					if( oe_C->len[2] + Nodes[A*2].len + oe_A->len[rev(i)] + Nodes[B*2-1].len + oe_B->len[rev(j)] <= Pairs[C].len )
						to_remove = true;
				}
			}
			else	// processing outArcs
			{
				if( oe_B->len[j] + Nodes[C*2-1].len + oe_C->len[2] + Nodes[A*2-1].len + oe_A->len[rev(i)] <= Pairs[B].len )
				{
					if( oe_C->len[2] + Nodes[A*2-1].len + oe_A->len[rev(i)] + Nodes[B*2-1].len + oe_B->len[rev(j)] <= Pairs[C].len )
						to_remove = true;
				}
			}
		}
	}
	else if(orient & OUTIE && orient2 & INNIE)
	{
		if( oe_A->len[i] + Nodes[B*2-1].len + oe_B->len[j] + Nodes[C*2].len + oe_C->len[2] <= Pairs[A].len)
		{
			if(d)	// processing inArcs
			{
				if( oe_B->len[j] + Nodes[C*2].len + oe_C->len[2] + Nodes[A*2].len + oe_A->len[rev(i)] <= Pairs[B].len )
				{
					if( oe_C->len[2] + Nodes[A*2].len + oe_A->len[rev(i)] + Nodes[B*2].len + oe_B->len[rev(j)] <= Pairs[C].len )
						to_remove = true;
				}
			}
			else	// processing outArcs
			{
				if( oe_B->len[j] + Nodes[C*2].len + oe_C->len[2] + Nodes[A*2-1].len + oe_A->len[rev(i)] <= Pairs[B].len )
				{
					if( oe_C->len[2] + Nodes[A*2-1].len + oe_A->len[rev(i)] + Nodes[B*2].len + oe_B->len[rev(j)] <= Pairs[C].len )
						to_remove = true;
				}
			}
		}
	}
	else
	{
		if( oe_A->len[i] + Nodes[B*2-1].len + oe_B->len[j] + Nodes[C*2-1].len + oe_C->len[2] <= Pairs[A].len)
		{
			if(d)	// processing inArcs
			{
				if( oe_B->len[j] + Nodes[C*2-1].len + oe_C->len[2] + Nodes[A*2].len + oe_A->len[rev(i)] <= Pairs[B].len )
				{
					if( oe_C->len[2] + Nodes[A*2].len + oe_A->len[rev(i)] + Nodes[B*2].len + oe_B->len[rev(j)] <= Pairs[C].len )
						to_remove = true;
				}
			}
			else	// processing outArcs
			{
				if( oe_B->len[j] + Nodes[C*2-1].len + oe_C->len[2] + Nodes[A*2-1].len + oe_A->len[rev(i)] <= Pairs[B].len )
				{
					if( oe_C->len[2] + Nodes[A*2-1].len + oe_A->len[rev(i)] + Nodes[B*2].len + oe_B->len[rev(j)] <= Pairs[C].len )
						to_remove = true;
				}
			}
		}
	}

	if( to_remove )
	{
		oe_C->flag(REMOVABLE);	// This edge is now removable, and the other two are supported
		oe_A->flag(SUPPORTED, i);
		oe_B->flag(SUPPORTED, rev(j));
	}
}

// Transitively reduce the edges
template <class Vertex, class Arc>
void reduce_graph(Vertex * nd_ptr, int num, ReadNode * Nodes)
{
	Arc * oe_A;
	Arc ** in_visited = new Arc*[num+1];
	int * in_flags = new int[num+1];
	Arc ** out_visited = new Arc*[num+1];
	int * out_flags = new int[num+1];

	for(int A=0; A<=num; A++)
	{
		in_visited[A] = NULL;
		out_visited[A] = NULL;
		in_flags[A] = -1;
		out_flags[A] = -1;
	}

	for(int A=1; A<=num; A++)
	{
		LinkedIter<Arc> iter_A(nd_ptr[A].outArcs);
		for(int d=0; d<2; d++)
		{
			if(d)	iter_A.change_list(nd_ptr[A].inArcs);

			// Step1: Mark the adjacent nodes as visited
			for(oe_A = iter_A.get_first(); oe_A != NULL; oe_A = iter_A.get_next())
			{
				int i = iter_A.get_flg();
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

			// Step2: Trace the out neighbours one by one and mark any third node that is transitively reached to be removed
			for(oe_A = iter_A.get_first(); oe_A != NULL; oe_A = iter_A.get_next())
			{
				int i = iter_A.get_flg();
				int B = oe_A->toRead[i];		// record the information of node B

				Direction orient = oe_A->arrow[i];

				LinkedIter<Arc> iter_B(nd_ptr[B].inArcs);
				if(orient & INNIE)	// if coming inside...
					iter_B.change_list(nd_ptr[B].outArcs);

				for(Arc * oe_B = iter_B.get_first(); oe_B != NULL; oe_B = iter_B.get_next())
				{
					int j = iter_B.get_flg();
					int C = oe_B->toRead[j];

					if( (in_visited[C] != NULL) && (oe_B->arrow[j] & INNIE) )
					{
						Arc * oe_C = in_visited[C];
						int k = in_flags[C];

						if(Nodes == NULL)
							check_read_reduction(oe_B, oe_A, oe_C, i, j, k);
						else
							check_mate_reduction(oe_B, oe_A, oe_C, i, j, d, Nodes, nd_ptr, orient, INNIE);


					}
					if( (out_visited[C] != NULL) && (oe_B->arrow[j] & OUTIE) )
					{
						Arc * oe_C = out_visited[C];
						int k = out_flags[C];

						if(Nodes == NULL)
							check_read_reduction(oe_B, oe_A, oe_C, i, j, k);
						else
							check_mate_reduction(oe_B, oe_A, oe_C, i, j, d, Nodes, nd_ptr, orient, OUTIE);
					}
					// Done with C, get next node
				}
				// Done with node B, so proceed with the next neighbour
			}

			// Step 3: Clean up
			for(oe_A = iter_A.get_first(); oe_A != NULL; oe_A = iter_A.get_next())
			{
				int i = iter_A.get_flg();
				in_visited[ oe_A->toRead[i] ] = NULL;
				in_flags[ oe_A->toRead[i] ] = -1;
				out_visited[ oe_A->toRead[i] ] = NULL;
				out_flags[ oe_A->toRead[i] ] = -1;
			}
		}
	}

	delete [] in_visited;
	delete [] out_visited;
	delete [] in_flags;
	delete [] out_flags;

	prune_graph<Vertex>(nd_ptr, num, REMOVABLE);
	LinkedNode<Arc>::shrink();
	Arc::shrink();
}

// removes all edges in the given graph (can be used to save space after chain collapsing)
template <class Vertex, class Arc>
void purge_edges(Vertex * nd_ptr, int num)
{
	if(nd_ptr != NULL)
	{
		for(int B=1; B<=num; B++)
			nd_ptr[B].remove_all();
	}
}

// prints stats about the graph (can be used for types inherited from BiNode)
template <class Vertex>
int check_graph(Vertex * nd_ptr, int num)
{
	int end_nodes = 0;
	int middlets = 0;
	int singlets = 0;
	int max_degree = 0;
	int degree;
	int num_nodes = 0;
	int num_arcs = 0;
	int others = 0;

	for(int x=1; x<=num; x++)
	{
		degree = nd_ptr[x].get_degree();
		num_arcs += degree;

		if(degree > max_degree)
			max_degree = degree;

		if(degree == 0)
			singlets++;
		else
		{
			num_nodes++;
			if(nd_ptr[x].is_end_node())
				end_nodes++;
			else if(nd_ptr[x].is_middle())
				middlets++;
			else
				others++;
		}
	}

	cout << "# Active nodes: " << num_nodes << endl;
	cout << "# End nodes: " << end_nodes << endl;
	cout << "# Middle nodes: " << middlets << endl;
	cout << "# Singleton nodes: " << singlets << endl;
	cout << "# Multi-degree nodes: " << others << endl;
	cout << "# arcs: " << num_arcs << endl;
	cout << "# max degree: " << max_degree << endl << endl;

	return num_nodes;	// number of active nodes
}

// can be used for all types inherited from BiNode
template <class Vertex>
int plot_distribution(char * filename, int num_of_blocks, int type, Vertex * nd_ptr, int num)
{
	ofstream fh;
	open_n_check(fh, filename);

	int A;
	int blocks[num_of_blocks+1];

	for(int i=0; i<=num_of_blocks; i++) {	blocks[i]=0; }

	for(A=1; A<=num; A++)
	{
		int deg = 0;

		switch(type)
		{
		case INDEG:
			deg = nd_ptr[A].get_in_degree();
			break;
		case OUTDEG:
			deg = nd_ptr[A].get_out_degree();
			break;
		case DEG:
			deg = nd_ptr[A].get_degree();
			break;
		default:
			throw "OvlGraph: Invalid type";
			break;
		}

		if(deg < num_of_blocks)
			blocks[deg]++;	// new kid on the block
		else
			blocks[num_of_blocks]++;

	}

	double x = 10.0;
	double y = 10.0;

	int coverage = 0;
	int highest = 0;

	fh << "graph [ " << "\n";
	for(int i=0; i<=num_of_blocks; i++)
	{
		if(blocks[i] > highest)
		{
			coverage = i;
			highest = blocks[i];
		}

		int width = ((blocks[i] * 1000) / num) + 2;
		fh << "node [ id " << i << " label " << i << "\n";
		fh << "LabelGraphics [ anchor \"w\" ] " << "\n";
		fh << "graphics [ x " << x - int(width/2) << " y " << y << " h 15.0 w " << width << "\n";
		fh << "fill \"#ffcc99\" line \"#000000\" type \"rectangle\" linewidth 1.0 ] ]" << "\n";
		y += 30;
	}
	check_n_close(fh);
	return coverage;
}

template <class Vertex, class Arc>
void chain_collapse(Vertex * nd_ptr, PathNode * new_ptr, int num)
{
	for(int A=1; A<=num; A++)
	{
		short int i, j;
		if( !(nd_ptr[A].is_middle()) )
		{
			LinkedIter<Arc> iter_A(nd_ptr[A].outArcs);
			Arc * e_ptr;

			for(int d=0; d<2; d++)
			{
				if(d)	iter_A.change_list(nd_ptr[A].inArcs);

				for(e_ptr = iter_A.get_first(); e_ptr != NULL; e_ptr = iter_A.get_next())
				{
					i = iter_A.get_flg();

					PathEdge * new_edge = new PathEdge();

					new_edge->toRead[0] = A;
					new_edge->arrow[0] = e_ptr->arrow[rev(i)];

					Arc * e_ptr2 = e_ptr;
					EdgeTuple * et;

					do
					{
						j = rev(i);

						et = new EdgeTuple(e_ptr2->toRead[i], e_ptr2->arrow[i], e_ptr2->len[i]);
						EdgeTuple * et2 = new EdgeTuple(e_ptr2->toRead[j], e_ptr2->arrow[j], e_ptr2->len[j]);

						new_edge->path[0].insert_entry(et2);
						new_edge->path[1].add_entry(et);

						new_edge->len += et2->offset;
						new_edge->len += et->offset;

						LinkedIter<Arc> iter( nd_ptr[ et->read ].inArcs );

						if(et->orientation & INNIE)
							iter.change_list( nd_ptr[ et->read ].outArcs );

						e_ptr2 = iter.get_first();
						if(e_ptr2 != NULL)
							i = iter.get_flg();

					} while( e_ptr2 != NULL && (nd_ptr[ et->read ].is_middle()) );

					new_edge->toRead[1] = et->read;
					new_edge->arrow[1] = (et->orientation & (INNIE|OUTIE));		// do not inherit other flags
					new_edge->len = new_edge->len / 2;

					// the latter is to make sure loops do not get written twice
					if(A < new_edge->toRead[1] || (A == new_edge->toRead[1] && new_edge->break_tie()) )
					{
						new_ptr[A].add_arc(new_edge, 1);
						new_ptr[ new_edge->toRead[1] ].add_arc(new_edge, 0);
					}
					else
						delete new_edge;
				}
			}
		}
	}
}

// prints the edges as contigs; use for a PathNode graph collapsed from ReadNodes only
void print_contigs(char * filename, PathNode * nd_ptr, int num, bool onestrand)
{
	std::ofstream graphFile;
	open_n_check(graphFile, filename);

	for(int A=1; A<=num; A++)
	{
		if(onestrand)
			nd_ptr[A].print_edges_ss(graphFile);
		else
			nd_ptr[A].print_edges(graphFile);
	}

	check_n_close(graphFile);
}

// removes nodes that have a single PathEdge with less than maxnum internal nodes
int remove_buds(PathNode * nd_ptr, int num, int maxnum)
{
	int thorns = 0;
	for(int x=1; x<=num; x++)
	{
		if( nd_ptr[x].get_degree() == 1 )
		{
			LinkedIter<PathEdge> iter(nd_ptr[x].inArcs);
			if(nd_ptr[x].get_out_degree() == 1)
				iter.change_list(nd_ptr[x].outArcs);

			PathEdge * pe = iter.get_first();
			int i = iter.get_flg();
			int y = pe->toRead[i];

			if(pe->path_size() <= maxnum)
			{
				if( ((pe->arrow[i] & INNIE) && nd_ptr[y].get_in_degree() > 1 && nd_ptr[y].get_out_degree() > 0) ||
					((pe->arrow[i] & OUTIE) && nd_ptr[y].get_out_degree() > 1 && nd_ptr[y].get_in_degree() > 0) )
				{
					pe->flag(REMOVABLE);
					nd_ptr[y].remove_marked(REMOVABLE);
					nd_ptr[x].remove_marked(REMOVABLE);
					thorns++;
				}
			}
		}
	}
	return thorns;
}

// --------------------------IN PROGRESS--------------------------------
// to do: add a check for conflicts and remove edges if necessary
int matesort(PathNode * nd_ptr, int x, int component, int * que, int * S, int * L)
{
	int qbegin = 0;
	int qend = 1;
	que[0] = x;
	nd_ptr[x].forward = true;
	nd_ptr[x].comp = component;

	while(qend > qbegin)
	{
		x = que[qbegin++];
		LinkedIter<PathEdge> iter(nd_ptr[x].inArcs);
		for(int d=0; d<2; d++)
		{
			if(d) iter.change_list(nd_ptr[x].outArcs);
			for(PathEdge * pe = iter.get_first(); pe != NULL; pe = iter.get_next())
			{
				int i = iter.get_flg();
				int y = pe->toRead[i];

				bool  before = nd_ptr[y].forward;

				if( (pe->arrow[0] & OUTIE) == (pe->arrow[1] & OUTIE) )		// they face the same direction
					nd_ptr[y].forward = !nd_ptr[x].forward;
				else
					nd_ptr[y].forward = nd_ptr[x].forward;

				if(nd_ptr[y].comp == component) // already in the list
				{
					if(before != nd_ptr[y].forward)
					{
						nd_ptr[y].forward = before;
						pe->flag(REMOVABLE);
						nd_ptr[y].remove_marked(REMOVABLE);
					}
					continue;
				}

				nd_ptr[y].comp = component;
				que[qend++] = y;
			}
		}
	}
	// now everything on the que has been marked with the component id and a fixed direction
	// next identify all the nodes that have no incoming edges
	int counter = 0;
	int sbegin = 0;
	int send = 0;
	for(int k=0; k<qend; k++)
	{
		int z = que[k];
		if( (nd_ptr[z].forward && nd_ptr[z].get_in_degree() == 0) || (!nd_ptr[z].forward && nd_ptr[z].get_out_degree() == 0) )
			S[send++] = z;
	}

	while(send > sbegin)
	{
		int z = S[sbegin++];
		L[counter++] = z;

		LinkedIter<PathEdge> iter(nd_ptr[z].outArcs);
		if(!nd_ptr[z].forward)
			iter.change_list(nd_ptr[z].inArcs);

		for(PathEdge * pe = iter.get_first(); pe != NULL; pe = iter.get_next())
		{
			int i = iter.get_flg();
			int y = pe->toRead[i];
			nd_ptr[y].count += 1;

			if( (nd_ptr[y].forward && nd_ptr[y].get_in_degree() == nd_ptr[y].count) ||
				(!nd_ptr[y].forward && nd_ptr[y].get_out_degree() == nd_ptr[y].count) )
			{
				S[send++] = y;
			}
		}
	}

	// now L is a topologically sorted list of the nodes in the component
	for(int k=0; k<counter; k++)
	{
		nd_ptr[ L[k] ].count = 0;	// reset so that we can use it as distance
		S[k] = 0;					// reset so that we can use it as the predecessor array
	}

	for(int k=0; k<counter; k++)
	{
		int z = L[k];

		LinkedIter<PathEdge> iter(nd_ptr[z].outArcs);
		if(!nd_ptr[z].forward)
			iter.change_list(nd_ptr[z].inArcs);

		for(PathEdge * pe = iter.get_first(); pe != NULL; pe = iter.get_next())
		{
			int i = iter.get_flg();
			int y = pe->toRead[i];

			if(nd_ptr[y].count <= nd_ptr[z].count + pe->path_size() )
			{
				nd_ptr[y].count = nd_ptr[z].count + pe->path_size();
				S[y] = z;	// update the predecessor of y
				que[y] = (pe->arrow[i] & (OUTIE|INNIE));
			}
		}
	}

	int maxi = 0;
	int winner = 0;
	for(int k=0; k<counter; k++)
	{
		if(nd_ptr[L[k]].count > maxi)
		{
			maxi = nd_ptr[L[k]].count;
			winner = L[k];
		}
	}
	return winner;		// this is the end node of the longest path
}

// removes bubbles of path size at most maxnum
int remove_bubbles(PathNode * nd_ptr, int num, int maxnum)
{
	int popped = 0;
	for(int x=1; x<=num; x++)
	{
		if(nd_ptr[x].get_degree() == 3 && !nd_ptr[x].is_end_node())
		{
			LinkedIter<PathEdge> iter(nd_ptr[x].inArcs);
			if(nd_ptr[x].get_out_degree() == 2)
				iter.change_list(nd_ptr[x].outArcs);

			PathEdge * pe = iter.get_first();
			int i = iter.get_flg();
			int y = pe->toRead[i];

			PathEdge * pe2 = iter.get_next();
			int j = iter.get_flg();
			int z = pe2->toRead[j];

			if( y==z && (pe->arrow[i] & (INNIE|OUTIE)) == (pe2->arrow[j] & (INNIE|OUTIE)) &&
				(pe->path_size() <= maxnum || pe2->path_size() <= maxnum))
			{
				PathEdge * pe3 = NULL;
				int k = -1;

				if(pe->path_size() < pe2->path_size() || pe->len < pe2->len)
				{
					pe3 = pe;
					k = i;
				}
				else
				{
					pe3 = pe2;
					k = j;
				}

				pe3->flag(REMOVABLE);

				nd_ptr[y].remove_marked(REMOVABLE);
				nd_ptr[x].remove_marked(REMOVABLE);
				popped++;
			}
		}
	}
	return popped;
}

#endif
