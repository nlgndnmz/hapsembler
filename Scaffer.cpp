/*****************************************************************************
    $Author: Nilgun Donmez$

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

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <climits>
#include <sstream>

#include "DnaRead.h"
#include "GraphDef.h"
#include "HapUtils.h"
#include "UGraph.h"
#include "ContigNode.h"
#include "SmithWaterman.h"

#include "lp_lib.h"			// required for lp_solve

#include "Scaffer.h"

using namespace std;

Scaffer::Scaffer(int maxDeg, int maxOvl, int minContigSize, int maxTransversal)
{
	numPairs = 0;
	numLibs = 0;
	numContigs = 0;

	maxOverlap = maxOvl;
	maxK = maxTransversal;
	minContig = minContigSize;
	maxDegree = maxDeg;

	Libs = NULL;
	Pairs = NULL;
}

Scaffer::~Scaffer()
{
	delete [] Libs;
	delete [] Pairs;
}

void Scaffer::read_info(char * filename, int minSupport)
{
	int st, end, insert, dev, ori;

	numLibs = 0;
	numPairs = 0;

	ifstream fh;
	open_n_check(fh, filename);
	while(fh >> st >> end >> insert >> dev >> ori)
	{
		numLibs++;
		if(end>numPairs)
			numPairs = end;
	}
	check_n_close(fh);

	numPairs /= 2;

	Pairs = new ReadPair[2*numPairs+1];
	Libs = new Library[numLibs+1];

	int counter = 0;

	ifstream fh2;
	open_n_check(fh2, filename);
	while(fh2 >> st >> end >> insert >> dev >> ori)
	{
		Libs[++counter].insertSize = insert;
		Libs[counter].deviation = (short int) (dev);
		Libs[counter].sdwidth = 3.0;
		Libs[counter].orient = 1;		// field reserved for future use

		Libs[counter].coverage = 0.0;
		Libs[counter].minSupport = minSupport;

		for(int i=st; i<=end; i++)
		{
			Pairs[i].lib = counter;
			Pairs[i].length = 1;
		}
	}
	check_n_close(fh2);
}

bool Scaffer::orient_component(int * que, int numMembers,
	int numEdges, int maxDegree, int * mapping, int comp)
{
	for(int k=0; k<numMembers; k++)
	{
		mapping[que[k]] = k;
		contigs[que[k]].component = comp;
	}

	numEdges *= 2;

	ContigEdge ** edgePointers = new ContigEdge*[numEdges];

	int dubnum = 2*numMembers;

	int numRemoved = 0;
	int numCut = 0;
	int edgeCounter = dubnum;
	int middleton = 100000;	// an edge that is the only incoming/outgoing edge for both nodes it connects, should rarely be removed
	int realNode = 200000;	// a real node, should only be removed after all edges are considered
	int aux2 = 300000;		// the second auxiliary node, need never be removed

	UGraph * ug = new UGraph(dubnum + numEdges, numEdges, maxDegree);
	for(int k=0; k<numMembers; k++)
	{
		int A = que[k];
		ug->add_node(2*mapping[A], realNode - contigs[A].get_degree());		// the more edges, the more likely a misassembly
		ug->add_node(2*mapping[A]+1, realNode - contigs[A].get_degree());
		ug->add_edge(2*mapping[A], 2*mapping[A]+1);

		LinkedIter<ContigEdge> iter(contigs[A].inArcs);
		int degree = contigs[A].get_in_degree();

		for(int d=0; d<2; d++)
		{
			if(d)
			{
				iter.change_list(contigs[A].outArcs);
				degree = contigs[A].get_out_degree();
			}

			for(ContigEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
			{
				int i = iter.get_flg();
				int C = e_ptr->toRead[i];

				if(A > C || contigs[A].component != comp || contigs[C].component != comp)
					continue;

				int source = (d) ? 2*mapping[A]+1 : 2*mapping[A];
				int target = (e_ptr->arrow[i] & OUTIE) ? 2*mapping[C]+1 : 2*mapping[C];

				int degree2 = (e_ptr->arrow[i] & INNIE) ? contigs[C].get_in_degree() : contigs[C].get_out_degree();

				int penalty = contigs[A].get_degree() + contigs[C].get_degree();

				if(degree2 == 1 && degree == 1)
					ug->add_node(edgeCounter, middleton + e_ptr->support - penalty);
				else
					ug->add_node(edgeCounter, e_ptr->support - penalty);

				ug->add_node(edgeCounter+1, aux2);

				ug->add_edge(source, edgeCounter);
				ug->add_edge(edgeCounter, edgeCounter+1);
				ug->add_edge(edgeCounter+1, target);

				edgePointers[edgeCounter-dubnum] = e_ptr;
				edgePointers[edgeCounter+1-dubnum] = e_ptr;
				edgeCounter+=2;
			}
		}
	}

	bool ok = true;
	int * X = new int[maxK+1];
	int nX = ug->oddcycle(X, maxK, ok);
	delete ug;

	if(ok)
	{
		for(int i=0; i<nX; i++)
		{
			if(X[i] < dubnum)
			{
				int A = (X[i]%2==0) ? X[i]/2 : (X[i]-1)/2;
				cout << "Removing contig " << que[A] << " (length: " << contigs[que[A]].len <<
					" contained: " << contigs[que[A]].containedLen << ")" << endl;
				contigs[que[A]].flag |= REMOVABLE;
				numRemoved++;
			}
			else	// it must be an edge
			{
				cout << "Removing the edge between " << edgePointers[X[i]-dubnum]->toRead[0] << " and "
					<< edgePointers[X[i]-dubnum]->toRead[1] << endl;
				edgePointers[X[i]-dubnum]->flag(REMOVABLE);
				numCut++;
			}
		}

		if(numCut > 0 || numRemoved > 0)
			cout << numCut << " edges and " << numRemoved << " contigs are found to be conflicting.." << endl;
	}

	delete [] X;
	delete [] edgePointers;

	return ok;
}


void Scaffer::mark_articulation_edges(int ** edges, ContigEdge ** edgePtrs, int numMembers, int n)
{
	int * stack = new int[n];
	int * depths = new int[n];
	int * parents = new int[n];
	int * lowpoints = new int[n];
	bool * inStack = new bool[n];

	for(int k=0; k<n; k++)
	{
		inStack[k] = false;
		lowpoints[k] = n+1;	// this serves as "infinity"
		depths[k] = n+1;
		parents[k] = -1;
	}

	int top = 0;
	int root = 0;

	stack[0] = root;
	depths[root] = top;
	inStack[root] = true;

	while(top >= 0)
	{
		int A = stack[top];
		bool found = false;
		bool cutVertex = false;
		int numKids = 0;
		int lowest = depths[A];

		for(int k=1; k<edges[A][0]; k++)
		{
			int C = edges[A][k];
			if(!inStack[C])
			{
				found = true;
				stack[++top] = C;		// push the stack
				depths[C] = top;
				inStack[C] = true;
				parents[C] = A;
				break;
			}
			else if(parents[A] != C)	// exclude the parent
				lowest = min(min(depths[C], lowpoints[C]), lowest);

			if(parents[C] == A)
			{
				if(A == root)
				{
					numKids++;
					if(numKids > 1)
						cutVertex = true;
				}
				else if(lowpoints[C] >= depths[A])
					cutVertex = true;
			}
		}

		if(!found)	// we have reached a dead end
		{
			top--;	// pop the stack
			lowpoints[A] = lowest;
			if(cutVertex == true && A > numMembers)
				edgePtrs[A-numMembers]->flag(NON_PROPER);
		}
	}

	delete [] stack;
	delete [] depths;
	delete [] parents;
	delete [] lowpoints;
	delete [] inStack;
}


void Scaffer::find_articulation_points(int * que, int * map, int numMembers)
{
	int maxDegree = 0;
	int numEdges = 0;
	for(int k=0; k<numMembers; k++)
	{
		int A = que[k];
		maxDegree = max( maxDegree, contigs[A].get_degree());
		numEdges += contigs[A].get_degree();
	}

	int n = numMembers + (numEdges/2);
	ContigEdge ** edgePtrs = new ContigEdge*[numEdges/2];
	int ** edges = new int*[n];
	for(int k=0; k<n; k++)
	{
		edges[k] = new int[maxDegree+1];
		edges[k][0] = 1;
	}

	int ind = numMembers;
	for(int k=0; k<numMembers; k++)
	{
		int A = que[k];
		LinkedIter<ContigEdge> iter(contigs[A].outArcs);
		for(int d=0; d<2; d++)
		{
			for(ContigEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
			{
				int i = iter.get_flg();
				int B = e_ptr->toRead[i];
				if(i == 1)
				{
					edges[ind][(edges[ind][0])++] = map[A];
					edges[ind][(edges[ind][0])++] = map[B];

					edges[map[A]][(edges[map[A]][0])++] = ind;
					edges[map[B]][(edges[map[B]][0])++] = ind;

					edgePtrs[ind-numMembers] = e_ptr;		// backdoor
					ind += 1;
				}
			}
			iter.change_list(contigs[A].inArcs);
		}
	}

	mark_articulation_edges(edges, edgePtrs, numMembers, n);

	for(int k=0; k<n; k++)
		delete [] edges[k];
	delete [] edges;
	delete [] edgePtrs;
}

// note: this routine assumes the graph is simple (i.e. no multi-edges)
bool Scaffer::process_biconnected_components(int * que, int numMembers, int * map, int * map2, int step, int lib, int & gain)
{
	if(step == 1 && numMembers < maxK+2)	// otherwise do not bother articulation
		return orient_component(que, numMembers, numMembers*numMembers, numMembers, map, 1);

	for(int k=0; k<numMembers; k++)
		map[que[k]] = k;

	find_articulation_points(que, map, numMembers);		// marks each cut edge as NON_PROPER

	bool * inStack = new bool[numMembers];
	int * parents = new int[numMembers];
	for(int k=0; k<numMembers; k++)
	{
		contigs[que[k]].component = -1;
		inStack[k] = false;
	}

	int comp = 1;
	bool completed = true;
	for(int k=0; k<numMembers; k++)
	{
		int A = que[k];
		if(inStack[map[A]])
			continue;

		parents[0] = A;
		inStack[map[A]] = true;

		int begin = 0;
		int end = 1;

		int numEdges = 0;
		int maxDegree = 0;

		while(end-begin > 0)
		{
			A = parents[begin++];
			numEdges += contigs[A].get_degree();
			if(contigs[A].get_degree() > maxDegree)
				maxDegree = contigs[A].get_degree();

			LinkedIter<ContigEdge> iter(contigs[A].outArcs);
			for(int d=0; d<2; d++)
			{
				for(ContigEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
				{
					int i = iter.get_flg();
					int B = e_ptr->toRead[i];

					if(!inStack[map[B]] && (!e_ptr->flagged(NON_PROPER)))
					{
						parents[end++] = B;
						inStack[map[B]] = true;
					}
				}
				iter.change_list(contigs[A].inArcs);
			}
		}

		if(end < 3)
			continue;

		if(step == 1)
		{
			if(!orient_component(parents, end, numEdges, maxDegree, map2, ++comp))
			{
				completed = false;
				break;
			}
		}
		else if(step == 2)
			feedback_arc_set(parents, end, ++comp);
		else if(step == 3)
		{
			int x = solve_lp(parents, end, map2, ++comp, lib);
			if(x < 0)
				cout << "lp_solve failed to return" << endl;
			else
				gain += x;
		}
	}

	delete [] inStack;
	delete [] parents;

	return completed;
}


void Scaffer::feedback_arc_set(int * que, int numMembers, int comp)
{
	int * s1 = new int[numMembers+1];
	int * s2 = new int[numMembers+1];
	int ind1 = 0;
	int ind2 = 0;

	for(int k=0; k<numMembers; k++)
		contigs[que[k]].component = comp;

	while(ind1 + ind2 < numMembers)
	{
		int maxDelta = -numMembers;
		int maxSupport = -1;
		int winner = 0;

		for(int t=0; t<3; t++)
		{
			int numRemoved = 1;
			while(numRemoved > 0)
			{
				numRemoved = 0;
				for(int k=0; k<numMembers; k++)
				{
					int A = que[k];
					if(contigs[A].component != comp)
						continue;

					int dPlus = 0;
					int dMinus = 0;
					int totalSupport = 0;
					LinkedIter<ContigEdge> iter(contigs[A].inArcs);
					for(int d=0; d<2; d++)
					{
						int degree = 0;
						int support = 0;
						for(ContigEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
						{
							int i = iter.get_flg();
							int C = e_ptr->toRead[i];
							if(contigs[C].component == comp)
							{
								degree++;
								support += e_ptr->support;
							}
						}
						if((d==1 && contigs[A].forward) || (d==0 && !contigs[A].forward))
						{
							dPlus = degree;
							totalSupport = support;
						}
						else
							dMinus = degree;
						iter.change_list(contigs[A].outArcs);
					}

					if(t==2)
					{
						if(dPlus - dMinus > maxDelta || (dPlus - dMinus == maxDelta && totalSupport > maxSupport))
						{
							winner = A;
							maxDelta = dPlus - dMinus;
							maxSupport = totalSupport;
						}
					}
					else if(t==0 && dPlus == 0)	// the node is a sink or singleton
					{
						s2[ind2++] = A;
						contigs[A].component = -1;
						numRemoved++;
					}
					else if(t==1 && dMinus == 0)
					{
						s1[ind1++] = A;
						contigs[A].component = -1;
						numRemoved++;
					}
				}
			}
		}

		if(winner > 0)
		{
			s1[ind1++] = winner;
			contigs[winner].component = -1;
		}
	}

	// combine the two sequences
	for(int k=ind2-1; k>=0; k--)
		s1[ind1++] = s2[k];

	// finally, remove the feedback edges
	for(int k=0; k<ind1; k++)
	{
		int A = s1[k];
		contigs[A].component = -2;

		LinkedIter<ContigEdge> iter(contigs[A].outArcs);
		if(!contigs[A].forward)
			iter.change_list(contigs[A].inArcs);

		for(ContigEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
		{
			int i = iter.get_flg();
			int C = e_ptr->toRead[i];

			if(contigs[C].component == -2)	// C is already assigned
				e_ptr->flag(REMOVABLE);
		}
	}

	delete [] s1;
	delete [] s2;
}

bool Scaffer::assign_components(int step, int & gain, int lib)
{
	int numRemoved = 0;
	for(int k=1; k<=numContigs; k++)
	{
		if(contigs[k].flag & REMOVABLE)
			contigs[k].remove(contigs);
		else
			numRemoved += contigs[k].remove_marked(REMOVABLE);
	}
	cout << "Removing " << numRemoved/2.0 << " edges.." << endl;

	int * que = new int[numContigs+1];
	int * mapping = new int[numContigs+1];
	int * mapping2 = new int[numContigs+1];
	bool * assigned = new bool[numContigs+1];

	for(int k=1; k<=numContigs; k++)
		assigned[k] = false;

	bool completed = true;
	for(int k=1; k<=numContigs; k++)
	{
		if(assigned[k])
			continue;

		que[0] = k;
		assigned[k] = true;
		contigs[k].forward = true;

		int begin = 0;
		int end = 1;

		while(end-begin > 0)
		{
			int A = que[begin++];
			contigs[A].fill_que(contigs, que, end, assigned);
		}

		if(end < 3)
			continue;

		if(!process_biconnected_components(que, end, mapping, mapping2, step, lib, gain))
		{
			completed = false;
			break;
		}
	}

	delete [] que;
	delete [] mapping;
	delete [] mapping2;
	delete [] assigned;

	return completed;
}

void Scaffer::adjust_libraries(ReadTag * tags, bool calibrate, long int totalBases,
	int NX, ContigTag * ctags, double epsilon)
{
	long int * libMean = new long int[numLibs+1];
	int * libSupport = new int[numLibs+1];

	int * success = new int[numLibs+1];
	int * failure = new int[numLibs+1];

	int histSize = 100000;
	int ** histogram = new int*[numLibs+1];

	for(int k=0; k<=numLibs; k++)
	{
		libMean[k] = 0;
		libSupport[k] = 0;
		success[k] = 0;
		failure[k] = 0;

		histogram[k] = new int[histSize+1];
		for(int d=0; d<=histSize; d++)
			histogram[k][d] = 0;
	}		

	for(int k=1; k<=numPairs; k++)
	{
		int j = 2*k;
		int i = j-1;
		if(tags[j].strand == INNIE)
		{
			i = 2*k;
			j = i-1;
		}

		int t = Pairs[i].lib;			

		if(tags[j].contigID >= 0 && tags[i].contigID >= 0)
			Libs[t].coverage += 1;

		if(tags[i].contigID == tags[j].contigID && tags[i].contigID >= 0 && ctags[tags[i].contigID].length > (NX - 2))
		{
			int distance = tags[j].offset - tags[i].offset;
			if(tags[i].strand != tags[j].strand && distance > 0)	// otherwise the orientation and/or mapping must be wrong
			{
				distance += Pairs[j].length;				
				if(distance < histSize)
				{
					histogram[t][distance] += 1;
					libSupport[t] += 1;
					success[t]++;
					libMean[t] += distance;
				}
				else
				{
					histogram[t][histSize] += 1;	// it is clearly an outlier
					failure[t]++;
				}
			}
			else
				failure[t]++;
		}			
	}
	
	cout << "processed the reads" << endl;

	for(int k=1; k<=numLibs; k++)
	{
		double confidence = (double)libSupport[k]/Libs[k].coverage;
		Libs[k].coverage = totalBases/Libs[k].coverage;	// this will actually serve as the arrival rate
		libMean[k] = int(libMean[k]/(double)libSupport[k]);

		cout << "--- Statistics for library #" << k;
		cout << " (support: " << libSupport[k] <<", confidence: "<< confidence << ") --- " << endl;
		if(NX < 2*libMean[k] || success[k] < 10)	// the estimate may not be reliable
		{
			cout << "WARNING: support is too low, library statistics will not be altered!" << endl;
			continue;
		}
			
		double ratio = success[k] / (double) (success[k] + failure[k]);
		Libs[k].minSupport = max(Libs[k].minSupport, int(log(epsilon) / log(1-ratio)));

		cout << "Success ratio: " << ratio << " (" << success[k] << "/" << failure[k]+success[k] << ")" << endl;
		cout << "Estimated coverage: " << (Libs[k].coverage/totalBases) * libMean[k] << endl;
		cout << "Minimum support is set to " << Libs[k].minSupport << endl;

		// what follows is a bootstrap method to adjust the library mean and deviation
		if(calibrate)
			bootstrap_library(libMean, libSupport, histogram, histSize, Libs, k);
	}

	for(int k=0; k<=numLibs; k++)
		delete [] histogram[k];

	delete [] histogram;
	delete [] libMean;
	delete [] libSupport;

	delete [] success;
	delete [] failure;
}

int Scaffer::fill_contigs(ReadTag * tags, int lib)
{
	int numAdded = 0;
	for(int k=1; k<=numPairs; k++)
	{
		int j = 2*k;
		int i = j-1;

		if(tags[i].contigID != tags[j].contigID && tags[i].contigID >= 0 && tags[j].contigID >= 0)
		{
			int l = Pairs[i].lib;
			if(l != lib)
				continue;

			int distance = 0;

			Direction jDir = INNIE;
			int jCon = tags[j].contigID;
			int jPar = contigs[jCon].parent;

			Direction iDir = INNIE;
			int iCon = tags[i].contigID;
			int iPar = contigs[iCon].parent;

			if(iPar == jPar)	// both contigs are on the same scaffold
				continue;

			if(lib>1 && ((iPar == iCon && contigs[iCon].len < Libs[lib-1].insertSize)
				|| (jPar == jCon && contigs[jCon].len < Libs[lib-1].insertSize)))
				continue;

			if(tags[j].strand == INNIE)
			{
				if(contigs[jCon].same)
				{
					jDir = OUTIE;
					distance += (contigs[jPar].containedLen - contigs[jCon].offset - tags[j].offset);
				}
				else
				{
					jDir = INNIE;
					distance += (contigs[jCon].len + contigs[jCon].offset - tags[j].offset);
				}
			}
			else if(tags[j].strand == OUTIE)
			{
				if(contigs[jCon].same)
				{
					jDir = INNIE;
					distance += (contigs[jCon].offset + tags[j].offset + Pairs[j].length);
				}
				else
				{
					jDir = OUTIE;
					distance += (contigs[jPar].containedLen - contigs[jCon].offset -
						contigs[jCon].len + tags[j].offset + Pairs[j].length);
				}
			}

			if(tags[i].strand == INNIE)
			{
				if(contigs[iCon].same)
				{
					iDir = OUTIE;
					distance += (contigs[iPar].containedLen - contigs[iCon].offset - tags[i].offset);
				}
				else
				{
					iDir = INNIE;
					distance += (contigs[iCon].len + contigs[iCon].offset - tags[i].offset);
				}
			}
			else if(tags[i].strand == OUTIE)
			{
				if(contigs[iCon].same)
				{
					iDir = INNIE;
					distance += (contigs[iCon].offset + tags[i].offset + Pairs[i].length);
				}
				else
				{
					iDir = OUTIE;
					distance += (contigs[iPar].containedLen - contigs[iCon].offset -
						contigs[iCon].len + tags[i].offset + Pairs[i].length);
				}
			}

			if((distance - maxOverlap) < (Libs[l].insertSize + Libs[l].sdwidth * Libs[l].deviation))
			{
				double d = 0.0;
				if(tags[i].weight < 0 && tags[j].weight < 0)
					d = 20.0;

				// double d = double(tags[i].weight + tags[j].weight) / (Pairs[i].length + Pairs[j].length);
				ContigEdge * edge = new ContigEdge(iPar, jPar, iDir, jDir, Libs[l].insertSize - distance, Libs[l].deviation, l, d);

				if( contigs[iPar].add_arc(edge, 1) )
				{
					contigs[jPar].add_arc(edge, 0);
					numAdded++;
				}
				else
					delete edge;	// edge not needed
			}
		}
	}
	return numAdded;
}

void Scaffer::write_scafftmp(char * outputFilename)
{
	string s(outputFilename);
	s += ".scafftmp";

	char * ofname = new char[s.size()+1];
	strcpy(ofname, s.c_str());

	ofstream fh;
	open_n_check(fh, ofname);

	for(int A=1; A<=numContigs; A++)
	{
		if(contigs[A].parent != A || contigs[A].contains.get_size() < 2)
			continue;

		LinkedIter<Containee> iter(contigs[A].contains);
		for(Containee * ce = iter.get_first(); ce != NULL; ce = iter.get_next())
		{
			int B = ce->contig;
			fh << ce->mean << " " << ((contigs[B].same) ? 2 : 1) << " " << B << " " << contigs[B].len << "\n";
		}
		fh << "-1 0" << endl;
	}

	check_n_close(fh);
	delete [] ofname;
}

void Scaffer::chain_collapse()
{
	bool secondRound = false;
	for(int d=0; d<2; d++)
	{
		if(d) secondRound = true;

		for(int A=1; A<=numContigs; A++)
		{
			if( contigs[A].parent != A || (!contigs[A].is_middle() && !contigs[A].is_single_end()))
				continue;

			bool isOut = true;
			if(contigs[A].is_single_end())
			{
				if(contigs[A].get_in_degree() == 1)
					isOut = false;
			}
			else
			{
				LinkedIter<ContigEdge> iter(contigs[A].outArcs);
				ContigEdge * e_ptr = iter.get_first();
				int B = e_ptr->toRead[ iter.get_flg() ];
				iter.change_list(contigs[A].inArcs);
				e_ptr = iter.get_first();
				int C = e_ptr->toRead[ iter.get_flg() ];

				bool isBhub = (!contigs[B].is_middle() && !contigs[B].is_single_end());
				bool isChub = (!contigs[C].is_middle() && !contigs[C].is_single_end());

				if(!secondRound && isBhub == isChub)
					continue;

				if(isBhub && !isChub)
					isOut = false;
			}

			// first include the node itself
			if(contigs[A].contains.get_size() == 0)	// if it's a singleton
			{
				Containee * ce = new Containee;
				ce->contig = A;
				ce->mean = 0;
				ce->variance = 1;
				ce->support = 1;

				contigs[A].contains.add_entry(ce);

				if(!isOut)
					contigs[A].same = false;
			}
			else
			{
				if(!isOut)	// then we have to reverse the order of the contigs, otherwise everything remains as it is
				{
					LinkedIter<Containee> iterA(contigs[A].contains);
					LinkedList<Containee> * tmpList = new LinkedList<Containee>(false);		// do not delete data, just delete pointers
					Containee * prv = iterA.get_first();
					Containee * ce = iterA.get_next();

					while(ce != NULL)
					{
						contigs[prv->contig].same = !(contigs[prv->contig].same);
						prv->mean = ce->mean;
						prv->variance = ce->variance;
						prv->support = ce->support;

						tmpList->insert_entry(prv);
						prv = ce;
						ce = iterA.get_next();
					}

					contigs[prv->contig].same = !(contigs[prv->contig].same);
					prv->mean = 0;
					prv->variance = 1;
					prv->support = 1;
					tmpList->insert_entry(prv);

					contigs[A].contains.clr_list(false);
					contigs[A].contains.append_list_shallow(tmpList);	// this makes a shallow copy
					delete tmpList;
				}
			}

			LinkedIter<ContigEdge> iter(contigs[A].outArcs);
			if(!isOut)
				iter.change_list(contigs[A].inArcs);

			int B = -1;
			ContigEdge * e_ptr = iter.get_first();
			while(e_ptr != NULL && !e_ptr->flagged(REMOVABLE))
			{
				int i = iter.get_flg();
				B = e_ptr->toRead[i];

				if(B == A || (!contigs[B].is_middle() && !contigs[B].is_single_end()))
					break;

				if(contigs[B].contains.get_size() == 0)	// if it's a singleton
				{
					Containee * ce = new Containee;
					ce->contig = e_ptr->toRead[i];
					ce->mean = e_ptr->len[0];
					ce->variance = e_ptr->len[1];
					ce->support = e_ptr->support;

					contigs[A].contains.add_entry(ce);

					if(e_ptr->arrow[i] & OUTIE)
						contigs[B].same = false;
				}
				else
				{
					if(e_ptr->arrow[i] & INNIE)
					{
						LinkedIter<Containee> iterB(contigs[B].contains);
						Containee * ce = iterB.get_first();
						ce->mean = e_ptr->len[0];
						ce->variance = e_ptr->len[1];
						ce->support = e_ptr->support;

						contigs[A].contains.append_list(&contigs[B].contains);	// perform deep copy (requires a copy constructor)
					}
					else
					{
						LinkedIter<Containee> iterB(contigs[B].contains);
						LinkedList<Containee> * tmpList = new LinkedList<Containee>(false);	// do not delete data, just delete pointers
						Containee * prv = iterB.get_first();
						Containee * ce = iterB.get_next();

						while(ce != NULL)
						{
							contigs[prv->contig].same = !(contigs[prv->contig].same);
							prv->mean = ce->mean;
							prv->variance = ce->variance;
							prv->support = ce->support;

							tmpList->insert_entry(prv);
							prv = ce;
							ce = iterB.get_next();
						}

						contigs[prv->contig].same = !(contigs[prv->contig].same);
						prv->mean = e_ptr->len[0];
						prv->variance = e_ptr->len[1];
						prv->support = e_ptr->support;
						tmpList->insert_entry(prv);

						contigs[A].contains.append_list(tmpList);
						delete tmpList;
					}
					contigs[B].contains.clr_list();
				}
				e_ptr->flag(REMOVABLE);

				if(e_ptr->arrow[i] & INNIE)
					iter.change_list(contigs[B].outArcs);
				else
					iter.change_list(contigs[B].inArcs);

				e_ptr = iter.get_first();
			}

			int pos = 0;
			LinkedIter<Containee> iter2(contigs[A].contains);
			for(Containee * ce = iter2.get_first(); ce != NULL; ce = iter2.get_next())
			{
				int B = ce->contig;
				pos += ce->mean;
				contigs[B].offset = pos;
				contigs[B].parent = A;
				pos += contigs[B].len;
			}
			contigs[A].containedLen = pos;
		}
	}
}

int Scaffer::peekaboo(bool strict, int minSupp)
{
	int counter = 0;
	for(int A=1; A<=numContigs; A++)
	{
		LinkedIter<ContigEdge> iter(contigs[A].outArcs);
		for(int d=0; d<2; d++)
		{
			if(d)
			{
				if(contigs[A].get_in_degree() != 2)
					continue;
				iter.change_list(contigs[A].inArcs);
			}
			else
			{
				if(contigs[A].get_out_degree() != 2)
					continue;
			}

			ContigEdge * oe1 = iter.get_first();
			int k1 = iter.get_flg();
			int B = oe1->toRead[k1];

			ContigEdge * oe2 = iter.get_next();
			int k2 = iter.get_flg();
			int C = oe2->toRead[k2];

			bool midB = ( ((oe1->arrow[k1] & INNIE) && contigs[B].get_in_degree() > 1) ||
						  ((oe1->arrow[k1] & OUTIE) && contigs[B].get_out_degree() > 1) );

			bool midC = ( ((oe2->arrow[k2] & INNIE) && contigs[C].get_in_degree() > 1) ||
						  ((oe2->arrow[k2] & OUTIE) && contigs[C].get_out_degree() > 1) );

			if((!strict || (midB && !midC)) && oe1->support < minSupp && oe2->support > 3*oe1->support)
			{
				oe1->flag(REMOVABLE);
				contigs[B].remove_marked(REMOVABLE);
				counter++;
			}
			else if((!strict || (!midB && midC)) && oe1->support > 3*oe2->support && oe2->support < minSupp)
			{
				oe2->flag(REMOVABLE);
				contigs[C].remove_marked(REMOVABLE);
				counter++;
			}
			contigs[A].remove_marked(REMOVABLE);
		}
	}
	return counter;
}

// for qsort (in non-decreasing order)
int comparePStruct(const void * a, const void * b)
{
	return ( (*(PStruct*)a).x[0] - (*(PStruct*)b).x[0] );
}

int ** Scaffer::floyd_warshall(int * que, int n, int * map, int comp)
{
	int ** paths = new int*[n];
	for(int i=0; i<n; i++)
	{
		paths[i] = new int[n];
		for(int j=0; j<n; j++)
			paths[i][j] = n+1;		// this serves as infinity because no path could be greater
	}

	for(int k=0; k<n; k++)
	{
		int A = que[k];
		LinkedIter<ContigEdge> iter(contigs[A].outArcs);
		if(!contigs[A].forward)
			iter.change_list(contigs[A].inArcs);

		paths[map[A]][map[A]] = 0;
		for(ContigEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
		{
			int i = iter.get_flg();
			int B = e_ptr->toRead[i];
			if(contigs[B].component == comp)
			{
				paths[map[A]][map[B]] = 1;
				paths[map[B]][map[A]] = 1;
			}
		}
	}

	for(int k=0; k<n; k++)
	{
		for(int i=0; i<n; i++)
		{
			for(int j=0; j<n; j++)
			{
				paths[i][j] = min(paths[i][j], paths[i][k] + paths[k][j]);
			}
		}
	}
	return paths;
}

// note: for this function, numEdges should be accurate and not just an upper bound
int Scaffer::solve_lp(int * que, int numMembers, int * map, int comp, int lib)
{
	if(numMembers > 30)		// too large to tackle
		return 0;

	for(int k=0; k<numMembers; k++)
	{
		map[que[k]] = k;
		contigs[que[k]].component = comp;
	}

	int maxDist = 0;					// first, we need to calculate some values
	int numEdges = 0;
	for(int k=0; k<numMembers; k++)
	{
		int A = que[k];
		maxDist += contigs[A].len;
		LinkedIter<ContigEdge> iter(contigs[A].outArcs);
		if(!contigs[A].forward)
			iter.change_list(contigs[A].inArcs);

		for(ContigEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
		{
			int i = iter.get_flg();
			int B = e_ptr->toRead[i];
			if(contigs[B].component == comp)
				numEdges++;
		}
	}
	maxDist *= 2;
	int gain = -1;

#ifdef LPSOLVER
	lprec * lp;
	int n = numMembers+numEdges;		// number of free variables
	lp = make_lp(0, n);
	if(lp == NULL)
		return -1;

	int * colNo = NULL;		// these are for the constraints
	REAL * row = NULL;
	colNo = (int *) malloc(4 * sizeof(*colNo));
	row = (REAL *) malloc(4 * sizeof(*row));

	int * c = NULL;			// these are for the objective function
	REAL * r = NULL;
	c = (int *) malloc((n+1) * sizeof(*c));
	r = (REAL *) malloc((n+1) * sizeof(*r));

	if(colNo == NULL || row == NULL || r == NULL || c == NULL)
		return -1;

	for(int i=numMembers+1; i<=n; i++)			// set upper and lower bounds for the slack variables
		set_bounds(lp, i, 0.0, 1.0);		// note that by default all variables have a lower bound of 0.0

	set_add_rowmode(lp, TRUE);

	int ind = numMembers+1;					// note that column and row indices start from 1 in lp_solve
	int m = 0;
	for(int k=0; k<numMembers; k++)
	{
		int A = que[k];
		LinkedIter<ContigEdge> iter(contigs[A].outArcs);
		if(!contigs[A].forward)
			iter.change_list(contigs[A].inArcs);

		for(ContigEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
		{
			int i = iter.get_flg();
			int B = e_ptr->toRead[i];

			if(contigs[B].component != comp)
				continue;

			int j = 0;					// 1*map[A] -1*map[B] +maxDist*I_ab <= -D_ab +maxDist
			colNo[j] = map[A]+1;
			row[j++] = 1;
			colNo[j] = map[B]+1;		// we add +1 because the columns start at 1 in lp_solve
			row[j++] = -1;
			colNo[j] = ind;
			row[j++] = maxDist;

			if(!add_constraintex(lp, j, row, colNo, LE, maxDist - (e_ptr->len[0] + contigs[A].len)))
				return -1;

			j = 0;						// 1*map[B] -1*map[A] +maxDist*I_ab <= D_ab +maxDist
			colNo[j] = map[B]+1;
			row[j++] = 1;
			colNo[j] = map[A]+1;
			row[j++] = -1;
			colNo[j] = ind;
			row[j++] = maxDist;

			if(!add_constraintex(lp, j, row, colNo, LE, maxDist + (e_ptr->len[0] + contigs[A].len)))
				return -1;

			c[m] = ind++;					// populate the objective function
			r[m++] = e_ptr->support;
		}
	}

	set_add_rowmode(lp, FALSE); 			// turn off the model building
    if(!set_obj_fnex(lp, m, r, c))		// set the model objective
		return -1;

	set_maxim(lp);							// maximize the objective
    set_verbose(lp, IMPORTANT);				// be concise

    int ret = solve(lp);
    if(ret == OPTIMAL)
    {
    	gain = 0;
		//cout << "Objective value: " << get_objective(lp) << endl;
		get_variables(lp, r);

		// first, flag all previous edges as removable
		for(int k=0; k<numMembers; k++)
		{
			int A = que[k];
			LinkedIter<ContigEdge> iter(contigs[A].outArcs);
			if(!contigs[A].forward)
				iter.change_list(contigs[A].inArcs);

			for(ContigEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
			{
				int i = iter.get_flg();
				int B = e_ptr->toRead[i];
				if(contigs[B].component == comp)
				{
					e_ptr->flag(REMOVABLE);
					gain++;
				}
			}
		}

		// next, sort the contigs based on their positions
		PStruct * pos = new PStruct[numMembers];
		for(int j=0; j<numMembers; j++)
		{
			pos[j].x[0] = int(r[j]);					// contig position
			pos[j].x[1] = que[j];					// contig ID
			pos[j].x[2] = contigs[que[j]].len;	// contig length
		}
		qsort(pos, numMembers, sizeof(PStruct), comparePStruct);

		int cutoff = 2;
		int ** paths = floyd_warshall(que, numMembers, map, comp);

		// last, add new edges based on the implicit order of the contigs
		for(int j=0; j<numMembers-1; j++)
		{
			int A = pos[j].x[1];
			int B = pos[j+1].x[1];

			int dist = pos[j+1].x[0] - (pos[j].x[0] + pos[j].x[2]);

			if(paths[map[A]][map[B]] > cutoff)	// they should not be adjacent
			{
				for(int k=j+2; k<numMembers; k++)
				{
					if(paths[map[A]][map[pos[k].x[1]]] <= cutoff)
					{
						B = pos[k].x[1];
						dist = pos[k].x[0] - (pos[j].x[0] + pos[j].x[2]);
						break;
					}
				}
			}

			Direction d1 = (contigs[A].forward) ? OUTIE : INNIE;
			Direction d2 = (contigs[B].forward) ? INNIE : OUTIE;

			ContigEdge * ce = new ContigEdge(A, B, d1, d2, dist, Libs[lib].deviation, lib, 0.0);
			ce->support = Libs[lib].minSupport + 1;		// we can not trust this edge too much
			if( contigs[A].add_fake_arc(ce, 1) )
				contigs[B].add_fake_arc(ce, 0);
			else
				delete ce;

			gain--;
		}
		delete [] pos;
		for(int j=0; j<numMembers; j++)
			delete [] paths[j];
		delete [] paths;
    }

	if(row != NULL)			// clean up
		free(row);
	if(colNo != NULL)
		free(colNo);
	if(r != NULL)
		free(r);
	if(c != NULL)
		free(c);
	if(lp != NULL)
		delete_lp(lp);
#endif

	return gain;
}

int Scaffer::prune_graph()
{
	int total = 0;
	for(int A=1; A<=numContigs; A++)
		total += contigs[A].remove_marked(REMOVABLE);		// remove all edges that are marked with a_flag
	return total;
}

long int Scaffer::read_mappings(char * inputFilename, char * outputFilename, bool calibrate)
{
	ifstream inputFile;
	open_n_check(inputFile, inputFilename);

	ContigTag * ctags = new ContigTag[2*numPairs+1];
	ReadTag * tags = new ReadTag[2*numPairs+1];

	for(int i=0; i<=2*numPairs; i++)
	{
		tags[i].contigID = -1;
		tags[i].strand = 0;
		tags[i].offset = 0;
		tags[i].weight = 0;

		ctags[i].length = 0;
		ctags[i].count = 0;
	}

	numContigs = 0;
	int read;
	while(inputFile >> read)
	{
		int contigNum, contigLength, offset, strand, len, weight;
		if(read == -1)	// then it means we have a new contig
		{
			inputFile >> contigNum >> contigLength;
			ctags[contigNum].length = contigLength;
			if(contigNum > numContigs)	numContigs = contigNum;
		}
		else
		{
			inputFile >> strand >> contigNum >> offset >> len >> weight;
			Pairs[read].length = len;

			if(Pairs[read].lib > 0 && ctags[contigNum].length > minContig)
			{
				if(tags[read].contigID == -1 || (tags[read].weight >= 0 && weight < 0))
				{
					tags[read].contigID = contigNum;
					tags[read].offset = offset;
					tags[read].weight = weight;

					if( Libs[Pairs[read].lib].orient == 3 || (Libs[Pairs[read].lib].orient == 2 && read%2==0) )
						tags[read].strand = (strand == 2) ? 1 : 2;
					else
						tags[read].strand = strand;
				}
				else	// then the read is not trustable
					tags[read].contigID = -2;
			}
			ctags[contigNum].count += Pairs[read].length;	// we use sequence coverage so we can combine multiple libraries
		}
	}
	check_n_close(inputFile);

	long int totalBases = 0;
	int * contigLens = new int[numContigs];
	for(int k=1; k<=numContigs; k++)
	{
		contigLens[k-1] = ctags[k].length;
		totalBases += ctags[k].length;
	}
	
	int N50 = length_stats(contigLens, numContigs, totalBases, 0.5);
	adjust_libraries(tags, calibrate, totalBases, N50, ctags, 0.001);		// re-calculate the library mean and variance

    contigs = new ContigNode[numContigs+1];
	for(int k=1; k<=numContigs; k++)
	{
		contigs[k].id = k;
		contigs[k].parent = k;
		contigs[k].len = ctags[k].length;
		contigs[k].containedLen = ctags[k].length;
		contigs[k].coverage = (ctags[k].count / (double)ctags[k].length);
	}
	
	cout << "starting the loop" << endl;

	for(int lib=1; lib<=numLibs; lib++)
	{
		maxOverlap = Libs[lib].insertSize;
		cout << "Number of initial contig edges: " << fill_contigs(tags, lib) << endl;

		for(int k=1; k<=numContigs; k++)
			contigs[k].finalize(Libs[lib].minSupport, contigs, Libs);		// marks weakly supported edges as removable
		cout << prune_graph() << " edges were removed due to low support" << endl;

		for(int k=1; k<=numContigs; k++)
			contigs[k].remove_conflicts(maxDegree);			// removes conflicting edges
		cout << prune_graph() << " edges were removed due to conflicts" << endl;

		int dummy = 0;
		while(!assign_components(1, dummy))
		{
			int maxDeg = 0;
			for(int k=1; k<=numContigs; k++)
			{
				int num = contigs[k].get_degree();
				if(num > maxDeg) maxDeg = num;
			}
			cout << "Decreasing maximum degree to " << maxDeg << endl;

			for(int k=1; k<=numContigs; k++)
				contigs[k].prune_nodes(contigs, maxDeg);
			prune_graph();
		}
		assign_components(2, dummy);	// this call assigns each contig as forward or reverse and removes cyclic edges
		cout << prune_graph() << " edges were removed due to cycles" << endl;
		cout << peekaboo(true, 2*Libs[lib].minSupport) << " edges were removed due to forks" << endl;

		int gain = 1;
		while(gain > 0)
		{
			gain = 0;
			assign_components(3, gain, lib);
			prune_graph();
			gain += peekaboo(true, 2*Libs[lib].minSupport);
			cout << "Graph smoothing : " << gain << endl;
		}

		chain_collapse();
		for(int k=1; k<=numContigs; k++)
			contigs[k].remove_all();
		
		ContigEdge::shrink();
		LinkedNode<ContigEdge>::shrink();	
	}		

	write_scafftmp(outputFilename);
	
	cout << "has finished writing the scaffold templates" << endl; 

	delete [] ctags;
	delete [] tags;
	delete [] contigs;
	delete [] contigLens;

	ContigEdge::purge();
	LinkedNode<ContigEdge>::purge();
	LinkedNode<Containee>::purge();

	return totalBases;
}

int Scaffer::write_scaffolds(char * contigFilename, char * outputFilename, long int totalLength)
{
	ofstream outputFile;
	open_n_check(outputFile, outputFilename);

	string s(outputFilename);
	s += ".scafftmp";

	char * ifname = new char[s.size()+1];
	strcpy(ifname, s.c_str());

	ifstream scaffoldFile;
	open_n_check(scaffoldFile, ifname);
	
	DnaRead * contigSeq = new DnaRead[numContigs+1];
	char * sequences = new char[totalLength+numContigs+1];
	
	bool * written = new bool[numContigs+1];
	for(int i=0; i<=numContigs; i++)
		written[i] = false;

	int buffsize = 32*1024*1024;
	char * buff = new char[buffsize];
	long int pos = 0;
	bool defline = true;
	int currentContig = 0;

	char * contigName = new char[1024];
	char * contigNum = new char[1024];
	int nameCounter = 0;
	int sequenceCounter = 0;

	ifstream fastaFile;
	open_n_check(fastaFile, contigFilename);

	while(true)
	{
		fastaFile.read(buff, buffsize);
		int c = fastaFile.gcount();
		for(int i=0; i<c; i++)
		{
			if(buff[i] == '\n')
			{
				if(defline)
				{
					int k = 0;

					for(int j=7; j<nameCounter; j++)
						contigNum[k++] = contigName[j];
					contigNum[k] = '\0';
					nameCounter = 0;

					contigSeq[currentContig].length = sequenceCounter;	// finalize the previous contig
					pos += sequenceCounter;

					currentContig = atoi(contigNum);
					contigSeq[currentContig].seq = &sequences[pos];			// assign the memory
					contigSeq[currentContig].seq[0] = '-';
					sequenceCounter = 1;
				}

				defline = false;
			}
			else if(buff[i] == '>')		// a new defline
				defline = true;
			else if(defline)
				contigName[nameCounter++] = buff[i];
			else
				contigSeq[currentContig].seq[sequenceCounter++] = buff[i];
		}
		if(!fastaFile) break;
	}

	contigSeq[currentContig].length = sequenceCounter;	// finalize the last contig
	check_n_close(fastaFile);
	
	int scaffoldNum = 0;
	int bytes = 0;
	int textwidth = 80;
	bool newScaffold = true;

	long int total = 0;
	int * scaffoldLengths = new int[numContigs+1];
	int counter = 0;

	DnaRead * rd1 = new DnaRead(maxOverlap+5);
	DnaRead * rd2 = new DnaRead(maxOverlap+5);

	SmithWaterman * SW = new SmithWaterman(maxOverlap+10, 10, 1, -2, -6, 0.06);

	int offset, orientation, contig, length;
	int prevContig = -1;

	int numOverlaps = 0;
	while(scaffoldFile >> offset >> orientation)
	{
		if(orientation == 0)	// a new scaffold
		{
			newScaffold = true;
			scaffoldLengths[counter++] = bytes;
			total += bytes;

			if(bytes%textwidth != 0) outputFile << "\n";
			bytes = 0;
		}
		else
		{
			if(newScaffold)
			{
				scaffoldNum += 1;
				outputFile << ">scaffold_" << scaffoldNum << "\n";
				newScaffold = false;
				contig = -1;
			}

			prevContig = contig;
			scaffoldFile >> contig >> length;

			written[contig] = true;

			if(orientation == OUTIE)
				contigSeq[contig].revcomp(true);
			else
				contigSeq[contig].revcomp(false);

			int st = 1;

			if(prevContig != -1 && offset < maxOverlap)	// then check whether there is an overlap
			{
				int m = maxOverlap;
				if(maxOverlap > contigSeq[prevContig].length)
					m = contigSeq[prevContig].length;

				if(maxOverlap > contigSeq[contig].length)
					m = contigSeq[contig].length;

				rd1->replicate_partial(contigSeq[prevContig], contigSeq[prevContig].length-m+1, contigSeq[prevContig].length);
				rd2->replicate_partial(contigSeq[contig], 1, m);

				bool isContained = false;
				SW->reset_cache();
				if( SW->align((*rd1), (*rd2), 1, 0-rd2->length, 2*rd2->length, false, isContained, false) )
				{
					OvlCache * oc = SW->get_cache();
					st = oc[0].pos[3] + 1;			// the other position should be zero anyways
					numOverlaps++;
					offset = 0;
				}
				else if(offset <= 0)
					offset = 10;	// just a fixed number of Ns to indicate there is a gap
			}

			for(int i=0; i<offset; i++)
			{
				outputFile << "N";
				bytes++;
				if(bytes%textwidth == 0) outputFile << "\n";
			}

			contigSeq[contig].write_seq(outputFile, st, contigSeq[contig].length-1, bytes, textwidth);
		}
	}
	check_n_close(scaffoldFile);

	cout << numOverlaps << " contigs had overlaps" << endl;

	for(int i=0; i<=numContigs; i++)
	{
		if(contigSeq[i].length > 0 && !written[i])
		{
			outputFile << ">contig_" << i << "\n";
			bytes = 0;
			contigSeq[i].write_seq(outputFile, 1, contigSeq[i].length-1, bytes, textwidth);
			scaffoldLengths[counter++] = bytes;
			total += bytes;

			if(bytes%textwidth != 0) outputFile << "\n";
		}
		contigSeq[i].detach();
	}
	check_n_close(outputFile);

	length_stats(scaffoldLengths, counter, total);

	delete [] buff;
	delete [] contigSeq;
	delete [] sequences;
	delete [] ifname;
	delete [] contigName;
	delete [] contigNum;
	delete [] written;
	delete [] scaffoldLengths;

	delete rd1;
	delete rd2;
	delete SW;

	return scaffoldNum;
}
