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

#include <fstream>
#include <iostream>

#include "EdgeTuple.h"
#include "GraphDef.h"

// define the static variables
EdgeTuple ** EdgeTuple::pool = 0;
int EdgeTuple::poolnum = 0;
int EdgeTuple::chunkctr = 0;

void * EdgeTuple::operator new(size_t sz)
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
		pool[++poolnum] = new EdgeTuple[chunksize];
		chunkctr = 1;
		return &((pool[poolnum])[chunkctr-1]);
	}
	else
		throw "EdgeTuple pool is out of memory!";
}

void EdgeTuple::operator delete(void * p)
{
	// don't do anything
}

void EdgeTuple::purge()
{
	std::cout << "EdgeTuple pool is being purged" << std::endl;
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

void EdgeTuple::init()
{
	std::cout << "EdgeTuple pool is being created" << std::endl;
	pool = new EdgeTuple*[poolsize];
	pool[0] = new EdgeTuple[chunksize];
	chunkctr = 0;
}

// generic constructor
EdgeTuple::EdgeTuple()
{
	read = 0;
	orientation = INNIE;
	offset = 0;
}

// full constructor
EdgeTuple::EdgeTuple(int read, Direction orientation, EdgeLength offset)
{
	this->read = read;
	this->orientation = orientation;
	this->offset = offset;
}

// copy constructor
EdgeTuple::EdgeTuple(const EdgeTuple & et)
{
	read = et.read;
	orientation = et.orientation;
	offset = et.offset;
}

EdgeTuple::~EdgeTuple()
{
	// nothing to do
}

void EdgeTuple::pickle(std::ofstream & fh)
{
	fh << offset << " " << orientation << " " << read << "\n";
}

void EdgeTuple::unpickle(std::ifstream & fh)
{
	fh >> offset >> orientation >> read;
}

