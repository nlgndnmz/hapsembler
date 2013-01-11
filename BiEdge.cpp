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

#include "BiEdge.h"

BiEdge::BiEdge()
{
	arrow[0] = NONE;
	arrow[1] = NONE;
	toRead[0] = 0;
	toRead[1] = 0;
}

BiEdge::~BiEdge()
{
	// nothing to delete
}

void BiEdge::flag(EdgeFlag f)
{
	arrow[0] |= f;
	arrow[1] |= f;
}

void BiEdge::flag(EdgeFlag f, int which)
{
	if(which!=0 && which!=1)
		throw "BiEdge: invalid index to arrow field";

	arrow[which] |= f;
}

bool BiEdge::flagged(EdgeFlag f)
{
	if((arrow[0] & f) && (arrow[1] & f))
		return true;
	return false;
}
