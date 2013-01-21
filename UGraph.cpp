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

#include <iostream>
#include <fstream>

#include "UGraph.h"
#include "DirNode.h"
#include "UNode.h"
#include "UEdge.h"
#include "HapUtils.h"

UGraph::UGraph(int num, int numE, int maxE)
{
	numNodes = num;
	nodes = new UNode[num];
	dirnodes = new DirNode[num*2+2];
	diredges = new DirEdge[16*(numE+num)];

	numEdges = 0;

	for(int i=0; i<num*2+2; i++)
		dirnodes[i].init(16*maxE);
}

UGraph::~UGraph()
{
	for(int A=0; A<numNodes; A++)
	{
		LinkedIter<UEdge> iter(nodes[A].edges);

		UEdge * ue = iter.get_first();
		while( ue != NULL )
		{
			int j = iter.get_flg();
			int B = ue->toNode[j];

			if(A > B)
				ue = iter.del_cur(true);
			else
				ue = iter.del_cur();
		}
	}

	delete [] nodes;
	delete [] diredges;
	delete [] dirnodes;

	LinkedNode<UEdge>::purge();
}

void UGraph::add_edge(int source, int target)
{
	UEdge * ue = new UEdge(source, target);
	nodes[source].add_edge(ue, 1);
	nodes[target].add_edge(ue, 0);
}

void UGraph::add_node(int source, int weight)
{
	dirnodes[source].weight = weight;
	dirnodes[source+numNodes].weight = weight;
}

void UGraph::add_diredge(int from, int to)
{
	diredges[numEdges].to = to;
	diredges[numEdges].from = from;
	diredges[numEdges].flow = 0;

	dirnodes[to].add_edge(diredges[numEdges]);
	dirnodes[from].add_edge(diredges[numEdges]);

	numEdges += 1;
}

void UGraph::augment_X_edges(int * X, int nX)
{
	for(int i=0; i<nX; i++)
	{
		int A = X[i];
		LinkedIter<UEdge> iter(nodes[A].edges);
		for(UEdge * ue = iter.get_first(); ue != NULL; ue = iter.get_next())
		{
			int j = iter.get_flg();
			int B = ue->toNode[j];

			if(nodes[B].status==removed)
			{
				ue->spin[0] = 1;
				ue->spin[1] = 2;
			}
			else if(nodes[B].status==part1)
			{
				ue->spin[j] = 0;
				ue->spin[(j+1)%2] = 1;
			}
			else if(nodes[B].status==part2)
			{
				ue->spin[j] = 0;
				ue->spin[(j+1)%2] = 2;
			}
		}
	}
}

void UGraph::augment_Y_edges(int * Y, int nY)
{
	for(int i=0; i<nY; i++)
	{
		int D = Y[i];
		LinkedIter<UEdge> iter(nodes[D].edges);

		if(nodes[D].status == subsetA)
		{
			add_diredge(2*numNodes, D);				// add an edge from the source to D
			add_diredge(D+numNodes, 2*numNodes+1);		// add an edge to the sink from D's sibling
		}
		else
		{
			add_diredge(2*numNodes, D+numNodes);		// add an edge from the source to D's sibling
			add_diredge(D, 2*numNodes+1);				// add an edge to the sink from D
		}

		for(UEdge * ue = iter.get_first(); ue != NULL; ue = iter.get_next())
		{
			int j = iter.get_flg();
			int C = ue->toNode[j];

			if(nodes[C].status==blank || nodes[C].status==removed)
				continue;

			if(nodes[D].status == subsetA)
			{
				if(nodes[C].status==part2)
				{
					add_diredge(D, C);
				}
				else if(nodes[C].status==part1)
				{
					add_diredge(C+numNodes, D+numNodes);
				}
				else if(nodes[C].status==subsetA)			// otherwise no edge is possible
				{
					if(ue->spin[j]==2)
						add_diredge(D, C+numNodes);
					if(ue->spin[j]==1)
						add_diredge(C, D+numNodes);
				}
			}
			else	// this time D_0 is in partB and D_1 is in partA
			{
				if(nodes[C].status==part2)
				{
					add_diredge(C+numNodes, D);
				}
				else if(nodes[C].status==part1)
				{
					add_diredge(D+numNodes, C);
				}
				else if(nodes[C].status==subsetB)			// otherwise no edge is possible
				{
					if(ue->spin[j]==2)
						add_diredge(C+numNodes, D);
					if(ue->spin[j]==1)
						add_diredge(D+numNodes, C);
				}
			}
		}
	}
}

void UGraph::augment_S_edges(int * S1, int nS1, int * S2, int nS2)
{
	for(int i=0; i<nS1; i++)
	{
		int A = S1[i];
		add_diredge(A, A+numNodes);		// these edges are to ensure vertex disjoint paths

		LinkedIter<UEdge> iter(nodes[A].edges);

		for(UEdge * ue = iter.get_first(); ue != NULL; ue = iter.get_next())
		{
			int j = iter.get_flg();
			int B = ue->toNode[j];

			if(nodes[B].status==part2)
			{
				add_diredge(A+numNodes, B);
				add_diredge(B+numNodes, A);
			}
		}
	}

	for(int i=0; i<nS2; i++)
	{
		int A = S2[i];
		add_diredge(A, A+numNodes);
	}
}

void UGraph::clear_dirnodes(int t)
{
	for(int i=0; i<=t; i++)				// no need to go further since those nodes haven't been explored yet
	{
		dirnodes[i].reset();
		dirnodes[i+numNodes].reset();
	}

	dirnodes[2*numNodes].reset();		// delete the source edges
	dirnodes[2*numNodes+1].reset();		// delete the sink edges
	numEdges = 0;
}

int UGraph::ford_fulkerson(int demand, int * que)
{
	int totalFlow = 0;
	int source = 2*numNodes;
	int sink = 2*numNodes+1;

	while(totalFlow < demand)	// it is not possible to send more flow than demand
	{
		que[0] = source;
		int begin = 0;
		int end = 1;
		dirnodes[source].enqued = true;

		while(end-begin>0)		// while the que is not empty
		{
			int A = que[begin];
			if(A == sink)		// stop since we have reached the sink
				break;
			begin += 1;
			dirnodes[A].fill_que(A, que, dirnodes, end);
		}

		for(int i=0; i<end; i++)
			dirnodes[que[i]].enqued = false;

		if(end-begin>0)		// means there is indeed a path from the source to the sink with available capacity
		{
			totalFlow += 1;			// the size of the flow always increases by one
			int B = sink;
			while(B != source)		// now visit each node on the path from the sink to the source and update flow
			{
				int parent = dirnodes[B].parent;
				int num = dirnodes[B].num;
				if( dirnodes[parent].edges[num]->to == B )
					dirnodes[parent].edges[num]->flow = 1;	// send flow
				else
					dirnodes[parent].edges[num]->flow = 0;	// push back the flow
				B = parent;
			}
		}
		else
			break;
	}
	return totalFlow;
}

int Dcompare(const void * a, const void * b)
{
	return ((DirWeight*)a)->weight - ((DirWeight*)b)->weight;
}

int UGraph::solve_cut(int demand, int * W, int * que, DirWeight * Q2, int t)
{
	int totalFlow = ford_fulkerson(demand, que);
	if(totalFlow == demand)
		return totalFlow;	// no need to fill W

	int begin = 0;
	int end = 1;
	int source = 2*numNodes;
	int sink = source + 1;
	Q2[0].name = source;
	dirnodes[source].enqued = true;

	while(end-begin>0)
	{
		int A = Q2[begin++].name;
		dirnodes[A].fill_reachable(A, Q2, dirnodes, end);
	}

	for(int i=0; i<end; i++)
		dirnodes[Q2[i].name].enqued = false;

	qsort(Q2, end, sizeof(DirWeight), Dcompare);

	int counter = 0;
	int flow = totalFlow;
	for(int i=0; i<end; i++)
	{
		if(Q2[i].name == sink || Q2[i].name == source)	// source and sink do not count
			continue;

		dirnodes[source].reset_flow();
		dirnodes[sink].reset_flow();
		for(int j=0; j<=t; j++)
		{
			dirnodes[j].reset_flow();
			dirnodes[j+numNodes].reset_flow();
		}

		dirnodes[Q2[i].name].forbidden = true;
		if(ford_fulkerson(demand, que) < flow)		// then this vertex must be in the cut
		{
			W[counter++] = (Q2[i].name > numNodes) ? Q2[i].name - numNodes : Q2[i].name;
			flow -= 1;
		}
		else
			dirnodes[Q2[i].name].forbidden = false;
	}
	assert(counter == totalFlow);
	return counter;
}

void UGraph::repartition(int * S1, int & nS1, int * S2, int & nS2, int * que)
{
	nS1 = 0;	// reset the lists
	nS2 = 0;

	for(int i=0; i<numNodes; i++)
	{
		if(nodes[i].status==part1 || nodes[i].status==part2)
			nodes[i].status = part;
	}

	for(int i=0; i<numNodes; i++)
	{
		if(nodes[i].status == part)		// not assigned to a partition yet
		{
			que[0] = i;
			nodes[i].status = part1;
			S1[nS1++] = i;

			int begin = 0;
			int end = 1;

			while(end-begin>0)
			{
				int A = que[begin++];
				LinkedIter<UEdge> iter(nodes[A].edges);
				for(UEdge * ue = iter.get_first(); ue != NULL; ue = iter.get_next())
				{
					int j = iter.get_flg();
					int B = ue->toNode[j];

					if(nodes[B].status == part)
					{
						if(nodes[A].status == part1)
						{
							nodes[B].status = part2;
							S2[nS2++] = B;
						}
						else
						{
							nodes[B].status = part1;
							S1[nS1++] = B;
						}
						que[end++] = B;
					}
					assert(nodes[A].status != nodes[B].status);		// sanity check
				}
			}
		}
	}
}

// Instead of recursion, we will simply build the graph by including one vertex at a time
int UGraph::oddcycle(int * X, int maxK, bool & ok)
{
	int k = 0;
	int * Y = new int[maxK+1];
	int * W = new int[maxK+1];
	int * temp = new int[maxK+1];

	for(int i=0; i<k; i++)			// discard the first k vertices
	{
		X[i] = i;
		nodes[i].status = removed;
	}

	int nX = k;

	int * S1 = new int[numNodes];
	int * S2 = new int[numNodes];

	S1[0] = k;
	S2[0] = k+1;

	int nS1 = 1;
	int nS2 = 1;

	nodes[k].status = part1;		// the partitions initially have one node each
	nodes[k+1].status = part2;

	int maxAll = 1 << (maxK+1);
	bool ** subsets = new bool*[maxAll];	// this will be helpful when enumerating all subsets of X
	for(int i=0; i<maxAll; i++)			// the empty set is excluded
	{
		subsets[i] = new bool[maxK+1];

		subsets[i][0] = (i%2==0) ? false : true;
		for(int j=1; j<maxK+1; j++)
			subsets[i][j] = (i & (1 << j)) ? true : false;
	}

	int * que = new int[2*numNodes+2];		// these are to help with ford-fulkerson
	DirWeight * Q2 = new DirWeight[2*numNodes+2];

	int all = 1 << (k+1);
	for(int t=k+2; t<numNodes; t++)
	{
		X[nX++] = t;		// add a new node to the set
		nodes[t].status = removed;

		if(nX <= k)			// continue until we reach the limit
			continue;

		augment_X_edges(X, nX);	// augment the edges based on the partitions

		bool ok2continue = false;

		for(int i=1; i<all; i++)	// for all non-empty subsets of X
		{
			int counter = 0;
			for(int j=0; j<k+1; j++)	// populate the subset Y
			{
				if(subsets[i][j])
					Y[counter++] = X[j];
			}

			if(counter < 1)		// skip if Y is empty
				continue;

			int total = 1 << counter;
			for(int n=0; n<total; n++)		// for all valid partitions of Y into yA and yB
			{
				for(int m=0; m<counter; m++)	// form the partition
					nodes[Y[m]].status = (subsets[n][m]) ? subsetA : subsetB;

				augment_Y_edges(Y, counter);
				augment_S_edges(S1, nS1, S2, nS2);

				int nW = solve_cut(counter, W, que, Q2, t);	// solve max flow with demand "counter"

				clear_dirnodes(t);

				if(nW < counter)	// Voila! We found a transversal that is smaller than X
				{
					for(int m=0; m<counter; m++)
						nodes[Y[m]].status = part1;		// this arbitrary assignment is only temporary

					int c = 0;
					for(int m=0; m<k+1; m++)
					{
						if(nodes[X[m]].status == removed)
							temp[c++] = X[m];
					}

					for(int m=0; m<nW; m++)
						temp[c++] = W[m];

					for(int m=0; m<c; m++)
					{
						X[m] = temp[m];
						nodes[X[m]].status = removed;
					}

					nX = c;
					repartition(S1, nS1, S2, nS2, que);     // now re-partition the bipartite graph into S1 and S2

					ok2continue = true;
					break;
				}
			}

			if(ok2continue)
				break;

			for(int m=0; m<counter; m++)
				nodes[Y[m]].status = removed;
		}

		if(!ok2continue)	// failed to find a transversal of size k
		{
			if(k == maxK)	// reached the tolerance just return with failure
			{
				ok = false;
				break;
			}
			t -= 1;		// take a step back
			k += 1;		// increase the traversal size
			all = 1 << (k+1);
		}
	}

	// clean up
	delete [] Y;
	delete [] W;
	delete [] temp;
	delete [] S1;
	delete [] S2;
	delete [] que;
	delete [] Q2;

	for(int i=0; i<maxAll; i++)
		delete [] subsets[i];
	delete [] subsets;

	return nX;
}
