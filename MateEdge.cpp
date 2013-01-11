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
#include "MateEdge.h"

// define the static variables
MateEdge ** MateEdge::pool = 0;
int MateEdge::poolnum = 0;
int MateEdge::chunkctr = 0;

void * MateEdge::operator new(size_t sz)
{
	if(pool == 0)
		init();

	if(chunkctr < chunksize)	// i.e. if there is still some space in the current chunk
	{
		chunkctr += 1;
		return &((pool[poolnum])[chunkctr-1]);
	}
	else if((poolnum+1) < poolsize)		// then allocate the next chunk
	{
		pool[++poolnum] = new MateEdge[chunksize];
		chunkctr = 1;
		return &((pool[poolnum])[chunkctr-1]);
	}
	else
		throw "MateEdge pool is out of memory!";
}

void MateEdge::operator delete(void * p)
{
	// don't do anything
}

void MateEdge::shrink()
{
	std::cout << "shrinking MateEdge pool" << std::endl;
	if(pool != 0)
	{
		for(int i=0; i<poolnum; i++)
		{
			bool alldead = true;
			for(int j=0; j<chunksize; j++)
			{
				if(pool[i][j].alive)
					alldead = false;
			}
			if(alldead)
			{
				delete [] pool[i];
				pool[i] = 0;
			}
		}
	}
}

void MateEdge::purge()
{
	std::cout << "MateEdge pool is being purged" << std::endl;
	if(pool != 0)
	{
		for(int i=0; i<=poolnum; i++)
			delete [] pool[i];

		delete [] pool;
	}
	pool = 0;
	poolnum = 0;
	chunkctr = 0;
}

void MateEdge::init()
{
	std::cout << "MateEdge pool is being created" << std::endl;
	pool = new MateEdge*[poolsize];
	pool[0] = new MateEdge[chunksize];
	chunkctr = 0;
}

// generic constructor
MateEdge::MateEdge()
{
	len[0] = 0;
	len[1] = 0;
	len[2] = 0;
	alive = true;
}

MateEdge::MateEdge(int rd1, int rd2, ArrowType a1, ArrowType a2, EdgeLength l1, EdgeLength l2, EdgeLength l3)
{
    toRead[0] = rd1;
    toRead[1] = rd2;

    arrow[0] = a1;
    arrow[1] = a2;

    len[0] = l3;	// this swap is intentional, makes inheriting easier!
    len[1] = l1;
    len[2] = l2;
    alive = true;
}

MateEdge::~MateEdge()
{
	alive = false;
}

// print the edge for contig
void MateEdge::print_edge(std::ofstream & fh, int k)
{
	fh << (int) len[k] << " " << (int) (arrow[k] & (INNIE|OUTIE)) << " " << toRead[k] << "\n";
}

// pickle the edge to the file
void MateEdge::pickle(std::ofstream & fh)
{
	fh << toRead[0] <<" "<< toRead[1] <<" "<< len[0] <<" "<< len[1] <<" "<< len[2] <<" "<<
	 (int) arrow[0] <<" "<< (int) arrow[1] <<"\n";
}

// unpickle the edge from the file
void MateEdge::unpickle(std::ifstream & fh)
{
	int a0, a1;
	fh >> toRead[0] >> toRead[1] >> len[0] >> len[1] >> len[2] >> a0 >> a1;

	arrow[0] = (ArrowType) a0;
	arrow[1] = (ArrowType) a1;
}
