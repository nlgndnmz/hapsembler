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
#include "OvlEdge.h"

// define the static variables
OvlEdge ** OvlEdge::pool = 0;
int OvlEdge::poolnum = 0;
int OvlEdge::chunkctr = 0;

// overloaded allocators
void * OvlEdge::operator new(size_t sz)
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
		pool[++poolnum] = new OvlEdge[chunksize];
		chunkctr = 1;
		return &((pool[poolnum])[chunkctr-1]);
	}
	else
		throw "OvlEdge pool is out of memory!";
}

void OvlEdge::operator delete(void * p)
{
	// don't do anything
}

void OvlEdge::shrink()
{
	std::cout << "shrinking OvlEdge pool" << std::endl;
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

void OvlEdge::purge()
{
	std::cout << "OvlEdge pool is being purged" << std::endl;
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

void OvlEdge::init()
{
	std::cout << "OvlEdge pool is being created" << std::endl;
	pool = new OvlEdge*[poolsize];
	pool[0] = new OvlEdge[chunksize];
	chunkctr = 0;
}

// generic constructor (the other members are intialized in the base class)
OvlEdge::OvlEdge()
{
	len[0] = 0;
	len[1] = 0;
	alive = true;
}

// full constructor
OvlEdge::OvlEdge(int read1, int read2, EdgeLength len1, EdgeLength len2, ArrowType arrow1, ArrowType arrow2)
{
	toRead[0] = read1;
	toRead[1] = read2;
	len[0] = len1;
	len[1] = len2;
	arrow[0] = arrow1;
	arrow[1] = arrow2;
	alive = true;
}

OvlEdge::~OvlEdge()
{
	alive = false;
}

// print the edge (used for contig reporting)
void OvlEdge::print_edge(std::ofstream & fh, int k)
{
	fh << len[k] << " " << (int) (arrow[k] & (INNIE|OUTIE)) << " " << toRead[k] << "\n";
}

// pickle the edge into a file for later retrieval
void OvlEdge::pickle(std::ofstream & fh)
{
	fh << toRead[0] <<" "<< toRead[1] <<" "<< len[0] <<" "<< len[1] <<" "<<
	 (int) arrow[0] <<" "<< (int) arrow[1] << "\n";
}

// unpickle the edge from the file
void OvlEdge::unpickle(std::ifstream & fh)
{
	int a0, a1;
	fh >> toRead[0] >> toRead[1] >> len[0] >> len[1] >> a0 >> a1;

	arrow[0] = (ArrowType) a0;
	arrow[1] = (ArrowType) a1;
}
