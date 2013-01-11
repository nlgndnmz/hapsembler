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

#include "ReadAln.h"

ReadAln::ReadAln()
{
	read0 = 0;
	read1 = 0;
	qual0 = 0;
	qual1 = 0;

	longest_gap = 0;
	weight = 0;
	start = 0;
	pos = 0;
	diff = 0;
}

ReadAln::~ReadAln()
{
	delete [] read0;
	delete [] read1;
	delete [] qual0;
	delete [] qual1;
}

void ReadAln::set(int size)
{
	read0 = new char[size];
	read1 = new char[size];
	qual0 = new short[size];
	qual1 = new short[size];
}
