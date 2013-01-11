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

#include <cmath>

#include "ScaffEdge.h"

ScaffEdge::ScaffEdge(int to)
{
	toContig = to;
	mean = 0;
	variance = 0;
}

ScaffEdge::ScaffEdge(ScaffEdge * src)
{
	this->toContig = src->toContig;
	this->mean = src->mean;
	this->variance = src->variance;
}


ScaffEdge::~ScaffEdge()
{
	// nothing to do
}

void ScaffEdge::add_fuzz(int m, int v)
{
	Fuzzy * f = new Fuzzy;
	f->mu = m;
	f->var = v;
	fuzzies.add_entry(f);
}

// to do: modify this to remove outliers or discard the entire edge if there are many outliers
void ScaffEdge::finalize()
{
	double p = 0.0;
	double q = 0.0;

	LinkedIter<Fuzzy> iter(fuzzies);
	for(Fuzzy * f = iter.get_first(); f != NULL; f = iter.get_next())
	{
		p += (double)f->mu / f->var;
		q += 1.0 / f->var;
	}
	mean = int(p / q);
	variance = int(pow(1.0/sqrt(q), 2));
}

