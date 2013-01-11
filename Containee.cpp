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

#include "Containee.h"
#include "GraphDef.h"

Containee::Containee()
{
	read = 0;
	same = true;
	begin = 0;
	end = 0;
}

Containee::~Containee()
{
	// nothing to do
}

Containee::Containee(int read, bool same, EdgeLength begin, EdgeLength end)
{
	this->read = read;
	this->same = same;
	this->begin = begin;
	this->end = end;
}

void Containee::pickle(std::ofstream & fh)
{
	fh << read <<" "<< (int)(same) <<" "<< (int)(begin) <<" "<< (int)(end) << "\n";
}

void Containee::unpickle(std::ifstream & fh)
{
	int r, s, b, e;
	fh >> r >> s >> b >> e;
	read = r;
	same = (s==0) ? false : true;
	begin = (EdgeLength) b;
	end = (EdgeLength) e;
}
