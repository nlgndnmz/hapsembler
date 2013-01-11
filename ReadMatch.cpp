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

#include "ReadMatch.h"

ReadMatch::ReadMatch(int m)
{
	num_hits = 0;
	max_hits = m;

	hits = new short int[max_hits];
	counts = new short int[max_hits];
}

ReadMatch::~ReadMatch()
{
	delete [] hits;
	delete [] counts;
}

void ReadMatch::reset()
{
	num_hits = 0;
}

void ReadMatch::compress(double epsilon, int q, int len1, int len2, bool nogap)
{
	epsilon = epsilon / 2.0;

	int cluster_num = 0;

	int next_cutoff = hits[0];
	int cluster_mean = hits[0];
	int total_hits = counts[0];
	int cluster_size = 1;

	int leeway = (cluster_mean<=0) ? (int)(((len2+cluster_mean)/2.0)*epsilon)+1 : (int)(((len1-cluster_mean)/2.0)*epsilon)+1;
	if(nogap) leeway = 3;

	for(int i=1; i<num_hits; i++)	// start from one because we've processed 0 already
	{
		if(hits[i] < (next_cutoff + leeway) && hits[i] > (next_cutoff - leeway))	// add it to the current cluster
		{
			cluster_mean += hits[i];
			total_hits += counts[i];
			cluster_size += 1;
		}
		else
		{
			cluster_mean = int(cluster_mean / (double)cluster_size);
			int ovl_size = ( cluster_mean <= 0 ) ? len2 + cluster_mean : len1 - cluster_mean;
			int tau =  int(ovl_size/(double)q) - int(2*ovl_size*epsilon);
			if(total_hits >= tau )
				hits[cluster_num++] = cluster_mean;

			next_cutoff = hits[i];
			cluster_mean = hits[i];
			total_hits = counts[i];
			cluster_size = 1;

			leeway = (cluster_mean<=0) ? (int)(((len2+cluster_mean)/2.0)*epsilon)+1 : (int)(((len1-cluster_mean)/2.0)*epsilon)+1;
			if(nogap) leeway = 3;
		}
	}

	cluster_mean = int(cluster_mean / (double)cluster_size);
	int ovl_size = ( cluster_mean <= 0 ) ? len2 + cluster_mean : len1 - cluster_mean;
	int tau =  int(ovl_size/(double)q) - int(2*ovl_size*epsilon);
	if(total_hits > tau)
		hits[cluster_num++] = cluster_mean;

	num_hits = cluster_num;
}

bool ReadMatch::get_next_pos(int & next, double epsilon, int & k1, int & k2, int len1, int len2, int & leeway, int minOverlap)
{
	while(next < num_hits)
	{
		int estOverlap = 0;
		int pos = hits[next++];
		if(pos <= 0)
		{
			leeway = (int)((len2 + pos) * epsilon) + 1;
			estOverlap = len2 + pos;

			k1 = 1;
			k2 = (0 - pos) - leeway;
		}
		else
		{
			leeway = (int)((len1 - pos) * epsilon) + 1;
			estOverlap = len1 - pos;

			k1 = pos - leeway;
			k2 = -2 * leeway;
			if(k1 <= 0)
			{
				k2 += (-k1);
				k1 = 1;
			}
		}
		leeway = 2*leeway+1;
		estOverlap += leeway;

		if(estOverlap >= minOverlap)
			return true;
	}
	return false;
}

bool ReadMatch::enter_hit(int index)
{
	bool first = false;
	if(num_hits == 0)	// first index
		first = true;

	bool skip = false;
	for(int i=0; i<num_hits; i++)
	{
		if(hits[i] == index)
		{
			counts[i] += 1;
			skip = true;
			break;
		}
	}

	if(!skip && num_hits < max_hits)	// not among the previously found hits and not reached the limit
	{
		hits[num_hits] = index;
		counts[num_hits++] = 1;
	}

	return first;
}
