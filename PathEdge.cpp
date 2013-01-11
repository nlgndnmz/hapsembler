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
#include "PathEdge.h"

// define the static variables
PathEdge ** PathEdge::pool = 0;
int PathEdge::poolnum = 0;
int PathEdge::chunkctr = 0;

void * PathEdge::operator new(size_t size)
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
		pool[++poolnum] = new PathEdge[chunksize];
		chunkctr = 1;
		return &((pool[poolnum])[chunkctr-1]);
	}
	else
		throw "PathEdge pool is out of memory!";
}

void PathEdge::operator delete(void * p)
{
	// don't do anything
}

void PathEdge::purge()
{
	std::cout << "PathEdge pool is being purged" << std::endl;
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

void PathEdge::init()
{
	std::cout << "PathEdge pool is being created" << std::endl;
	pool = new PathEdge*[poolsize];
	pool[0] = new PathEdge[chunksize];
	chunkctr = 0;
}

PathEdge::PathEdge()
{
	len = 0;
}

PathEdge::PathEdge(int read1, int read2, ArrowType arrow1, ArrowType arrow2, PathLength len)
{
    toRead[0] = read1;
    toRead[1] = read2;

    arrow[0] = arrow1;
    arrow[1] = arrow2;

    this->len = len;
}

PathEdge::~PathEdge()
{
	// nothing to do
}

// return the path size (note that both path[0] and path[1] should have the same size)
int PathEdge::path_size()
{
	return path[0].get_size();
}

// this function is used to break ties when the PathEdge is actually a self-loop (i.e. the edge starts and ends on the same node)
bool PathEdge::break_tie()
{
	EdgeTuple * et0 = path[0].get_root()->data;
	EdgeTuple * et1 = path[1].get_root()->data;

	if(et0->read < et1->read)
		return true;
	else if(et0->read > et1->read)
		return false;

	return (et0->orientation & INNIE);		// then the directions must be different
}

// this function is used to concatenate two edges. It performs deep copy so the old versions should be deleted by the caller
void PathEdge::add_path(LinkedList<EdgeTuple> & a_path, int k)
{
	LinkedIter<EdgeTuple> iter(a_path);

	for(EdgeTuple * et = iter.get_first(); et != NULL; et = iter.get_next())
	{
		EdgeTuple * et2 = new EdgeTuple(*et);
		path[k].add_entry(et2);
	}
}

// this function is used to copy a path except for the first tuple. Performs deep copy.
void PathEdge::copy_ex_first(LinkedList<EdgeTuple> & a_path, int k)
{
	LinkedIter<EdgeTuple> iter(a_path);
	EdgeTuple * et = iter.get_first();
	et = iter.get_next();

	while( et != NULL)
	{
		EdgeTuple * et2 = new EdgeTuple(*et);
		path[k].add_entry(et2);
		et = iter.get_next();
	}
}

// this function is used to copy a path except for the last tuple. Performs deep copy.
void PathEdge::copy_ex_last(LinkedList<EdgeTuple> & a_path, int k)
{
	LinkedIter<EdgeTuple> iter(a_path);
	EdgeTuple * et_prev = iter.get_first();
	EdgeTuple * et_next = iter.get_next();

	while(et_next != NULL)
	{
		EdgeTuple * et2 = new EdgeTuple(*et_prev);
		path[k].add_entry(et2);
		et_prev = et_next;
		et_next = iter.get_next();
	}
}

// print the path of the edge
void PathEdge::print_edge(std::ofstream & fh, int k)
{
	LinkedIter<EdgeTuple> iter(path[k]);
	for(EdgeTuple * et = iter.get_first(); et != NULL; et = iter.get_next())
		et->pickle(fh);
}

// print the path of the edge if there are at least "min" number of unvisited nodes in the path
bool PathEdge::print_edge(std::ofstream & fh, int k, int * reported, int min)
{
	int notdone = 0;
	bool first = false;
	bool second = false;
	bool report = true;

	LinkedIter<EdgeTuple> iter(path[k]);
	for(EdgeTuple * et = iter.get_first(); et != NULL; et = iter.get_next())
	{
		if(reported[et->read] == 0)
		{
			first = true;	// have seen our first unreported nodes
			notdone++;
		}
		if(first && reported[et->read] > 0)
			second = true;

		if(second && reported[et->read] == 0)
			report = false;		// there are breaks in the edge so do not report
	}

	if(report && notdone > min)
	{
		for(EdgeTuple * et = iter.get_first(); et != NULL; et = iter.get_next())
		{
			if(reported[et->read] == 0)
			{
				et->pickle(fh);
				reported[et->read] += 1;
			}
		}
		return true;
	}
	return false;
}

// pickle the edge into the file given
void PathEdge::pickle(std::ofstream & fh)
{
	fh<<(int)data_id<<" "<<path_size()<<" "<<toRead[0]<<" "<<toRead[1]<<" "<<(int)arrow[0]<<" "<<(int)arrow[1]<<" "<<len<<"\n";
	for(int d=0; d<2; d++)
		print_edge(fh, d);
}

// unpickle the edge from the file
void PathEdge::unpickle(std::ifstream & fh)
{
	int dummy, num, a0, a1;
	fh >> dummy >> num >> toRead[0] >> toRead[1] >> a0 >> a1 >> len;

	arrow[0] = (ArrowType) a0;
	arrow[1] = (ArrowType) a1;

	for(int d=0; d<2; d++)
	{
		for(int i=0; i<num; i++)
		{
			EdgeTuple * et = new EdgeTuple();
			et->unpickle(fh);
			path[d].add_entry(et);
		}
	}
}
