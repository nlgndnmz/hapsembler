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

#include <cstdlib>
#include <cmath>
#include <iostream>

#include "HapUtils.h"

using namespace std;

// opens an input file (ifstream)
void open_n_check(ifstream & fh, char * filename)
{
	fh.open(filename);
	if(!fh)
	{
		cerr << "HapUtils: Unable to open file (ifstream) " << filename << endl;
		exit(1);
	}
}

// opens an output file (ofstream)
void open_n_check(ofstream & fh, char * filename)
{
	fh.open(filename);
	if(!fh)
	{
		cerr << "HapUtils: Unable to open file (ofstream) " << filename << endl;
		exit(1);
	}
}

// closes an input file (ifstream)
void check_n_close(ifstream & fh)
{
	if(!fh.eof() && fh.fail())		// eof check is necessary since it also sets the fail flag
	{
		cerr << "HapUtils: file read failed" << endl;
		exit(1);
	}
	fh.close();
}

// closes an output file (ofstream)
void check_n_close(ofstream & fh)
{
	if(fh.fail())
	{
		cerr << "HapUtils: file write failed" << endl;
		exit(1);
	}
	fh.close();
}

// swap elements at indices qs1 and qs2 in que (also update parent pointers ptrs)
void swap_que_elts(int qs1, int qs2, QStruct * que, PStruct * ptrs)
{
	int readID = que[qs1].readID;
	Direction dir = que[qs1].dir;
	EdgeLength val = que[qs1].val;

	que[qs1].readID = que[qs2].readID;
	que[qs1].dir = que[qs2].dir;
	que[qs1].val = que[qs2].val;

	que[qs2].readID = readID;
	que[qs2].dir = dir;
	que[qs2].val = val;

	ptrs[ que[qs1].readID ].x[ que[qs1].dir ] = qs1;
	ptrs[ que[qs2].readID ].x[ que[qs2].dir ] = qs2;

}

// Required by qsort at various places, sorts in non-decreasing order
int compare(const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

int length_stats(int * lengths, int num, long int total, double frac)
{
	qsort(lengths, num, sizeof(int), compare);
	int nX = 0;
	long int soFar = 0;	
	long int cutoff = (int)(frac * total);
	for(int i=num-1; i>=0; i--)
	{
		soFar += lengths[i];
		if(soFar > cutoff)
		{
			cerr << "N" << int(frac*100) << ": " << lengths[i] << endl;
			nX = lengths[i];
			break;
		}
	}

	cerr << "Largest size: " << lengths[num-1] << endl;
	cerr << "Total number: " << num << endl;
	cerr << "Total length: " << total << endl;

	return nX;
}

// Required by the lexdfs method in OvlGraph.h
bool larger(const vector<int> & v1, const vector<int> & v2)	// is v1 lexicographically larger than v2?
{
	vector<int>::const_reverse_iterator it1;
	vector<int>::const_reverse_iterator it2;
    for( it1 = v1.rbegin(), it2 = v2.rbegin(); it1 != v1.rend() && it2 != v2.rend(); it1++, it2++ )
    {
    	if( (*it1) > (*it2) )	// then v1 is larger
			return true;
		else if( (*it1) < (*it2) )	// then v2 is larger
			return false;
    }
    // if we haven't returned yet, then either one the array is exhausted or they are equal
    if( it1 == v1.rend() )	// then v2 is larger (or equal)
		return false;

	return true;
}

// Required by the lexdfs method in OvlGraph.h
bool larger_eq(const vector<int> & v1, const vector<int> & v2)	// is v1 lexicographically larger than or equal to v2 ?
{
	vector<int>::const_reverse_iterator it1;
	vector<int>::const_reverse_iterator it2;
    for( it1 = v1.rbegin(), it2 = v2.rbegin(); it1 != v1.rend() && it2 != v2.rend(); it1++, it2++ )
    {
    	if( (*it1) > (*it2) )	// then v1 is larger
			return true;
		else if( (*it1) < (*it2) )	// then v2 is larger
			return false;
    }
    // if we haven't returned yet, then either one of the arrays is exhausted or they are equal
    if( it2 == v2.rend() )	// then v1 is larger or equal to v2
		return true;

	return false;
}

// used by correct_reads in Dna.h
void calculate_het_cutoff(int * table, int max_overlaps, double err, int max_t, double prob)
{
	table[0] = 1;
	table[1] = 2;
	for(int n=2; n<=max_overlaps; n++)
	{
		table[n] = max_t;
		for(int t=3; t<=max_t; t++)
		{
			double x = pow(err, t) * pow(1-err, n-t);
			long long int y = 1;
			double z = 1.0;

			for(int j=0; j<t; j++)
			{
				y = y * (n-j);
				z = z * (j+1);
			}
			if(x*(y/z) < prob)
			{
				table[n] = t;
				break;
			}
		}
	}
}

void bootstrap_library(long long int * libMean, int * libSupport, int ** histogram,
	int histSize, Library * Libs, int k, bool upperOnly)
{
	long long int newMean = 0;
	long long int newDev = 0;
	int ceiling = 0;
	double sdWidth = 0.0;

	for(int d=0; d<histSize; d++)
		newDev += (histogram[k][d] * (d-libMean[k]) * (d-libMean[k]));
	newDev = int(sqrt(newDev / (double)libSupport[k]));

	int total = 0;
	for(int d=0; d<histSize; d++)
	{
		total += histogram[k][d];
		if(total / (double)libSupport[k] > 0.99)
		{
			sdWidth = (d-libMean[k])/(double)newDev;
			ceiling = d;
			break;
		}
	}

	int iter = 1;
	while(iter < 100)
	{
		int n = libMean[k] - int(sdWidth*newDev);
		int m = (histSize < libMean[k] + int(sdWidth*newDev)) ? histSize : libMean[k] + int(sdWidth*newDev);

		newMean = 0;
		total = 0;
		for(int d=n; d<m; d++)					// calculate the mean using the bounded area only
		{										// if data is symmetrical around the mean this should not alter
			newMean += (d*histogram[k][d]);		// the mean too much
			total += histogram[k][d];
		}
		newMean = newMean/total;

		if(newMean - libMean[k] < 0.01*libMean[k] && libMean[k] - newMean < 0.01*libMean[k])
			break;

		libMean[k] = newMean;

		newDev = 0;
		for(int d=0; d<histSize; d++)
			newDev += (histogram[k][d] * (d-libMean[k]) * (d-libMean[k]));
		newDev = int(sqrt(newDev / (double)libSupport[k]));

		sdWidth = (ceiling - libMean[k])/(double)newDev;
		iter++;
	}

	double LB = 0.0;
	double UB = 0.0;

	total = 0;
	for(int d=0; d<histSize; d++)
	{
		total += histogram[k][d];
		if(total / (double)libSupport[k] > 0.99)
		{
			UB = (d-libMean[k])/(double)newDev;
			break;
		}
	}

	total = 0;
	for(int d=histSize-1; d>0; d--)
	{
		total += histogram[k][d];
		if(total / (double)libSupport[k] > 0.99)
		{
			LB = (libMean[k]-d)/(double)newDev;
			break;
		}
	}

	if(upperOnly)
		sdWidth = UB;
	else
		sdWidth = max(3.0, max(LB, UB));	// do not allow anything smaller

	cerr << "Adjusting mean insert size from " << Libs[k].insertSize << " to " << libMean[k] << endl;
	cerr << "Adjusting insert size deviation from " << Libs[k].deviation << " to " << newDev << endl;
	cerr << "Adjusting maximum standard deviations from " << Libs[k].sdwidth << " to " << sdWidth << endl;
	cerr << "( " << iter << " iterations were performed )" << endl;

	Libs[k].insertSize = (int) libMean[k];
	Libs[k].deviation = (short int) newDev;
	Libs[k].sdwidth = sdWidth;
}

