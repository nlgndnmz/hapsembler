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

#include "UNode.h"
#include "UEdge.h"

UNode::UNode()
{
	edges.set_no_delete();
	status = blank;
}

UNode::~UNode()
{
	// nothing to do
}

void UNode::add_edge(UEdge * ue, int which)
{
	edges.add_entry(ue, which);
}


void UNode::write2file(std::ofstream & fh, int limit)
{
	LinkedIter<UEdge> iter(edges);
	for(UEdge * ue = iter.get_first(); ue != NULL; ue = iter.get_next())
	{
		int j = iter.get_flg();
		int B = ue->toNode[j];
		if(B < limit)
			fh << ue->toNode[0] << " " << ue->toNode[1] << "\n";
	}
}
