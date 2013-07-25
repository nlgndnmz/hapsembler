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

#ifndef DNA_TEMPL_H
#define DNA_TEMPL_H

#include "LinkedIter.h"
#include "PathNode.h"
#include "GraphDef.h"

// removes node A from the graph
template <class Vertex, class Arc>
void remove_node(Vertex * nd_ptr, int A)
{
	LinkedIter<Arc> iter(nd_ptr[A].outArcs);
	for(int d=0; d<2; d++)
	{
		if(d)	iter.change_list(nd_ptr[A].inArcs);
		for(Arc * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
		{
			int i = iter.get_flg();
			int B = e_ptr->toRead[i];
			if(A != B)
				nd_ptr[B].remove_neighbour(A);
		}
	}
	nd_ptr[A].clear_edges();
}

// can be used for ReadNode and MatePair
template <class Vertex, class Arc>
int remove_buds(Vertex * nd_ptr, int num)
{
	int counter = 0;
	for(int A=1; A<=num; A++)
	{
		bool safe = false;
		if( nd_ptr[A].is_end_node() )	// if it is an end node
		{
			// this is a candidate for deleting, check if any neighbour has both directions
			LinkedIter<Arc> me_iter(nd_ptr[A].outArcs);
			if(nd_ptr[A].get_out_degree() == 0)
				me_iter.change_list(nd_ptr[A].inArcs);

			for(Arc * me = me_iter.get_first(); me != NULL; me = me_iter.get_next())
			{
				int i = me_iter.get_flg();
				int B = me->toRead[i];
				if((me->arrow[i] & OUTIE && nd_ptr[B].get_out_degree() > 1) || (me->arrow[i] & INNIE && nd_ptr[B].get_in_degree() > 1))
				{
					safe = true;	// then safe to delete A
					break;
				}
			}
		}

		if(safe)
		{
			remove_node<Vertex, Arc>(nd_ptr, A);
			counter++;
		}
	}
	return counter;
}

// use for ReadNode or MatePair only
template <class Vertex, class Arc>
int remove_contained_nodes(Vertex * nd_ptr, int num)
{
	int removed_nodes = 0;
	for(int A=1; A<=num; A++)
	{
		if(nd_ptr[A].flag & CONTAINED)
		{
			remove_node<Vertex, Arc>(nd_ptr, A);
			removed_nodes++;
		}
	}
	return removed_nodes;
}

// should be used for ReadNode and MatePair (use the overloaded function for PathNode)
template <class Vertex, class Arc>
int remove_bubbles(Vertex * nd_ptr, int num_vertices)
{
	int counter = 0;

	Arc ** in_visited = new Arc*[num_vertices+1];
	int * in_flags = new int[num_vertices+1];
	Arc ** out_visited = new Arc*[num_vertices+1];
	int * out_flags = new int[num_vertices+1];
	bool * winner = new bool[num_vertices+1];

	for(int x=0; x<=num_vertices; x++)
	{
		in_visited[x] = NULL;
		out_visited[x] = NULL;
		in_flags[x] = -1;
		out_flags[x] = -1;
		winner[x] = false;
	}

	int * ngbrs = new int[num_vertices+1];

	for(int x=1; x<num_vertices; x++)
	{
		if(nd_ptr[x].flag & DEAD)
			continue;

		int num = nd_ptr[x].get_neighbours(INNIE, ngbrs, false);
		for(int d=0; d<2; d++)
		{
			if(d) num = nd_ptr[x].get_neighbours(OUTIE, ngbrs, false);

			for(int i=0; i<num; i++)
			{
				bool remove = false;

				int A = abs(ngbrs[i]);
				if(nd_ptr[A].flag & DEAD)
					continue;

				LinkedIter<Arc> iter_A(nd_ptr[A].outArcs);
				for(int e=0; e<2; e++)
				{
					int neighbours = 0;
					if(e)	iter_A.change_list(nd_ptr[A].inArcs);

					Arc * oe_A;
					for(oe_A = iter_A.get_first(); oe_A != NULL; oe_A = iter_A.get_next())
					{
						int k = iter_A.get_flg();
						neighbours++;
						if(oe_A->arrow[k] & INNIE)
						{
							in_visited[ oe_A->toRead[k] ] = oe_A;
							in_flags[ oe_A->toRead[k] ] = k;
						}
						else
						{
							out_visited[ oe_A->toRead[k] ] = oe_A;
							out_flags[ oe_A->toRead[k] ] = k;
						}
					}

					for(int j=0; j<num; j++)
					{
						int B = abs(ngbrs[j]);
						if(A==B || nd_ptr[B].flag & DEAD)
							continue;

						bool same = (ngbrs[i] > 0 && ngbrs[j] >0) || (ngbrs[i] < 0 && ngbrs[j] < 0);
						int match = 0;

						LinkedIter<Arc> iter_B(nd_ptr[B].outArcs);
						if((e==1 && same) || (e==0 && !same))
							iter_B.change_list(nd_ptr[B].inArcs);

						for(Arc * oe_B = iter_B.get_first(); oe_B != NULL; oe_B = iter_B.get_next())
						{
							int k = iter_B.get_flg();
							int C = oe_B->toRead[k];

							if( (in_visited[C] != NULL) && (oe_B->arrow[k] & INNIE) )
								match++;
							if( (out_visited[C] != NULL) && (oe_B->arrow[k] & OUTIE) )
								match++;
						}

						if(match == neighbours)
						{
							if(e && winner[B])
							{
								remove = true;
								break;
							}
							else
								winner[B] = true;
						}
					}

					for(oe_A = iter_A.get_first(); oe_A != NULL; oe_A = iter_A.get_next())
					{
						int k = iter_A.get_flg();
						int B = oe_A->toRead[k];
						in_visited[ B ] = NULL;
						in_flags[ B ] = -1;
						out_visited[ B ] = NULL;
						out_flags[ B ] = -1;
					}

					if(e)
					{
						for(int j=0; j<num; j++)
						{
							int B = abs(ngbrs[j]);
							winner[B] = false;
						}
					}
				}

				if(remove)
				{
					remove_node<Vertex, Arc>(nd_ptr, A);
					counter++;
				}
			}
		}
	}

	delete [] ngbrs;
	delete [] winner;
	delete [] in_visited;
	delete [] out_visited;
	delete [] in_flags;
	delete [] out_flags;

	return counter;
}

int peekaboo(PathNode * paths, int num)
{
	int counter = 0;
	for(int A=1; A<=num; A++)
	{
		LinkedIter<PathEdge> iter_A(paths[A].outArcs);
		for(int d=0; d<2; d++)
		{
			if(d)
			{
				if(paths[A].get_in_degree() != 2)
					continue;
				iter_A.change_list(paths[A].inArcs);
			}
			else
			{
				if(paths[A].get_out_degree() != 2)
					continue;
			}

			PathEdge * oe1 = iter_A.get_first();
			int k1 = iter_A.get_flg();
			PathEdge * oe2 = iter_A.get_next();
			int k2 = iter_A.get_flg();

			if(oe1->path_size() > oe2->path_size() && oe2->path_size() == 1)
			{
				if( ((oe2->arrow[k2] & INNIE) && paths[oe2->toRead[k2]].get_in_degree() > 1) ||
					((oe2->arrow[k2] & OUTIE) && paths[oe2->toRead[k2]].get_out_degree() > 1) )
				{
					oe2->flag(REMOVABLE);
					paths[oe2->toRead[k2]].remove_marked(REMOVABLE);
					counter++;
				}
			}
			else if(oe2->path_size() > oe1->path_size() && oe1->path_size() == 1)
			{
				if( ((oe1->arrow[k1] & INNIE) && paths[oe1->toRead[k1]].get_in_degree() > 1) ||
					((oe1->arrow[k1] & OUTIE) && paths[oe1->toRead[k1]].get_out_degree() > 1) )
				{
					oe1->flag(REMOVABLE);
					paths[oe1->toRead[k1]].remove_marked(REMOVABLE);
					counter++;
				}
			}
			paths[A].remove_marked(REMOVABLE);
		}
	}
	return counter;
}

// generic method to syncronize bidirected edges/nodes with the collapsed version after simplifications
template <class Vertex, class Arc>
void path_sync(Vertex * nd_ptr, int num, PathNode * paths)
{
	bool * neighbours = new bool[num+1];
	for(int x=1; x<=num; x++)
	{
		neighbours[x] = false;
		LinkedIter<PathEdge> iter(paths[x].outArcs);
		for(int d=0; d<2; d++)
		{
			if(d)	iter.change_list(paths[x].inArcs);
			for(PathEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
			{
				int i = iter.get_flg();
				LinkedIter<EdgeTuple> iter(e_ptr->path[i]);
				for(EdgeTuple * et = iter.get_first(); et != NULL; et = iter.get_next())
					nd_ptr[et->read].flag |= VISITED;
			}
		}
	}

	for(int x=1; x<=num; x++)	// this loop will remove all the nodes absent from paths
	{
		if(paths[x].get_degree() == 0 && !(nd_ptr[x].flag & VISITED))
			remove_node<Vertex, Arc>(nd_ptr, x);
	}

	for(int x=1; x<=num; x++)	// this loop will remove all the edges absent from paths
	{
		if(paths[x].get_degree() == 0)	// which means it may be a middle node
			continue;

		LinkedIter<PathEdge> iter(paths[x].outArcs);
		LinkedIter<Arc> iter2(nd_ptr[x].outArcs);

		for(int d=0; d<2; d++)
		{
			if(d)
			{
				iter.change_list(paths[x].inArcs);
				iter2.change_list(nd_ptr[x].inArcs);
			}
			for(PathEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
			{
				int i = iter.get_flg();
				LinkedIter<EdgeTuple> iter3(e_ptr->path[i]);
				EdgeTuple * et = iter3.get_first();
				neighbours[et->read] = true;
			}
			for(Arc * e_ptr2 = iter2.get_first(); e_ptr2 != NULL; e_ptr2 = iter2.get_next())
			{
				int j = iter2.get_flg();
				if( !neighbours[e_ptr2->toRead[j]] )
					e_ptr2->flag(REMOVABLE);
			}
			for(PathEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
			{
				int i = iter.get_flg();
				LinkedIter<EdgeTuple> iter3(e_ptr->path[i]);
				EdgeTuple * et = iter3.get_first();
				neighbours[et->read] = false;
			}
		}
	}
	delete [] neighbours;
	prune_graph(nd_ptr, num, REMOVABLE);
}

// generic method to remove all edges from the graph matching the EdgeFlag a_flag
template <class Vertex>
long int prune_graph(Vertex * nd_ptr, int num, EdgeFlag a_flag)
{
	long int total = 0;
	for(int A=1; A<=num; A++)
		total += nd_ptr[A].remove_marked(a_flag);		// remove all edges that are marked with a_flag
	return total;
}

template <class Vertex, class Arc>
void reduce_graph(Vertex * nd_ptr, int num)
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
				int B = oe_A->toRead[i];

				Direction orient = oe_A->arrow[i];

				LinkedIter<Arc> iter_B(nd_ptr[B].inArcs);
				if(orient & INNIE)
					iter_B.change_list(nd_ptr[B].outArcs);

				for(Arc * oe_B = iter_B.get_first(); oe_B != NULL; oe_B = iter_B.get_next())
				{
					int j = iter_B.get_flg();
					int C = oe_B->toRead[j];

					Arc * oe_C = NULL;

					if( (in_visited[C] != NULL) && (oe_B->arrow[j] & INNIE) )
						oe_C = in_visited[C];

					if( (out_visited[C] != NULL) && (oe_B->arrow[j] & OUTIE) )
						oe_C = out_visited[C];

					if(oe_C != NULL)
					{
						int diff = oe_C->len[0] - (oe_B->len[0] + oe_A->len[0] + nd_ptr[B].len);
						if( abs(diff) < 6*(sqrt(oe_C->len[1]) + sqrt(oe_A->len[1]) + sqrt(oe_B->len[1])) )
						{
							oe_C->flag(REMOVABLE);		// flag as to be removed
							oe_A->flag(SUPPORTED, i);
							oe_B->flag(SUPPORTED, rev(j));
						}
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
}

// merges edges that become linear after various simplification procedures
void path_merge(PathNode * nd_ptr, int num)
{
	for(int x=1; x<=num; x++)
	{
		if( (nd_ptr[x].get_in_degree() == 1 && nd_ptr[x].get_out_degree() == 1) )
		{
			LinkedIter<PathEdge> iter(nd_ptr[x].inArcs);
			LinkedIter<PathEdge> iter2(nd_ptr[x].outArcs);

			PathEdge * pe = iter.get_first();
			int i = iter.get_flg();
			int y = pe->toRead[i];

			PathEdge * pe2= iter2.get_first();
			int j = iter2.get_flg();
			int z = pe2->toRead[j];

			if(y==x || z==x)	// do not merge loops
				continue;

			PathEdge * newpe = new PathEdge(y, z, pe->arrow[i], pe2->arrow[j], pe->len + pe2->len);

			newpe->add_path(pe2->path[rev(j)], 0);
			newpe->add_path(pe->path[i], 0);

			newpe->add_path(pe->path[rev(i)], 1);
			newpe->add_path(pe2->path[j], 1);

			nd_ptr[y].add_arc(newpe, 1);
			nd_ptr[z].add_arc(newpe, 0);

			pe2->flag(REMOVABLE);
			nd_ptr[z].remove_marked(REMOVABLE);
			nd_ptr[x].remove_marked(REMOVABLE);

			pe->flag(REMOVABLE);
			nd_ptr[y].remove_marked(REMOVABLE);
			nd_ptr[x].remove_marked(REMOVABLE);
		}
	}
}

#endif
