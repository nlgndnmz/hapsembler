/*****************************************************************************
    $Author: Nil $
    $Date: 2011-12-05 23:49:46 -0500 (Mon, 05 Dec 2011) $

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

#ifndef BI_EDGE_H
#define BI_EDGE_H

#include "GraphDef.h"

// A base class for a bidirected edge
class BiEdge
{
	public:

		BiEdge();
		~BiEdge();

		void flag(EdgeFlag);
		void flag(EdgeFlag, int);
		bool flagged(EdgeFlag);

		int toRead[2];
		ArrowType arrow[2];

	protected:

		const static int chunksize = 1048576;	// 1024x1024
		const static int poolsize = 1048576;	// altogether this will entail a terrabyte of items
};

#endif // BI_EDGE_H
