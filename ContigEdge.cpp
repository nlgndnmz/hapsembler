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

#include "ContigEdge.h"

// define the static variables
ContigEdge ** ContigEdge::pool = 0;
int ContigEdge::poolnum = 0;
int ContigEdge::chunkctr = 0;

void * ContigEdge::operator new(size_t sz)
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
		pool[++poolnum] = new ContigEdge[chunksize];
		chunkctr = 1;
		return &((pool[poolnum])[chunkctr-1]);
	}
	else
		throw "ContigEdge pool is out of memory!";
}

void ContigEdge::operator delete(void * p)
{
	if(chunkctr > 1 && p == &((pool[poolnum])[chunkctr-1]))	// if we are deleting the last added item, just reclaim the pointer
		chunkctr--;
}

void ContigEdge::shrink()
{
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

void ContigEdge::purge()
{
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

void ContigEdge::init()
{
	pool = new ContigEdge*[poolsize];
	pool[0] = new ContigEdge[chunksize];
	chunkctr = 0;
}

ContigEdge::ContigEdge()
{
	len[0] = 0;
	len[1] = 0;
	support = 0;
	lib = 0;
	weight = 0;
}

ContigEdge::ContigEdge(int contig0, int contig1, ArrowType arrow0, ArrowType arrow1, int mean, int sd, int libNum, double w)
{
	toRead[0] = contig0;
	toRead[1] = contig1;

	arrow[0] = arrow0;
	arrow[1] = arrow1;

	len[0] = mean;			// we will abuse these fields to use for mean and variance
	len[1] = sd * sd;
	support = 1;
	lib = libNum;
	weight = w;
}

ContigEdge::ContigEdge(const ContigEdge & ce)
{
	toRead[0] = ce.toRead[0];
	toRead[1] = ce.toRead[1];

	arrow[0] = ce.arrow[0];
	arrow[1] = ce.arrow[1];

	len[0] = ce.len[0];
	len[1] = ce.len[1];
	support = ce.support;
	lib = ce.lib;
	weight = ce.weight;
}

ContigEdge::~ContigEdge()
{
	// nothing to do
}

void ContigEdge::merge(ContigEdge * ce)
{
	this->support += 1;
	this->len[0] += ce->len[0];
	this->weight += ce->weight;
}

void ContigEdge::finalize(int minSupport)
{
	len[0] = int(len[0] / support);
	len[1] = int(len[1] / support);
	weight /= support;

	if(support < minSupport)
		this->flag(REMOVABLE);
}

void ContigEdge::print_edge(std::ofstream & fh)
{
	fh << toRead[0] <<" "<< toRead[1] <<" "<< ((arrow[0] & INNIE) ? "i" : "o") << ((arrow[1] & INNIE) ? "i " : "o ")
		<< (int) len[0] <<" "<< (int) len[1] <<" "<< weight <<" "<< support << "\n";
}
