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

#include "ScaffEdge.h"
#include "ScaffNode.h"

ScaffNode::ScaffNode()
{
	id = 0;
	len = 0;
	forward = true;

	parent = 0;
	offset = 0;
	offsetVar = 0;

	totalLen = 0;
	totalVar = 0;

	contained.set_no_delete();	// important!!
}

ScaffNode::~ScaffNode()
{
	// nothing to do
}

int ScaffNode::add_in_edge(int B, int mu, int sd)
{
	LinkedIter<ScaffEdge> iter(inArcs);
	for(ScaffEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
	{
		if(B == e_ptr->toContig)
		{
			e_ptr->add_fuzz(mu, sd);
			return 0;
		}
	}
	ScaffEdge * se = new ScaffEdge(B);
	se->add_fuzz(mu, sd);
	inArcs.add_entry(se);
	return 1;
}

int ScaffNode::add_out_edge(int B, int mu, int sd)
{
	LinkedIter<ScaffEdge> iter(outArcs);
	for(ScaffEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
	{
		if(B == e_ptr->toContig)
		{
			e_ptr->add_fuzz(mu, sd);
			return 0;
		}
	}
	ScaffEdge * se = new ScaffEdge(B);
	se->add_fuzz(mu, sd);
	outArcs.add_entry(se);
	return 1;
}

void ScaffNode::finalize()
{
	LinkedIter<ScaffEdge> iter(outArcs);
	for(int d=0; d<2; d++)
	{
		if(d) { iter.change_list(inArcs); }
		for(ScaffEdge * e_ptr = iter.get_first(); e_ptr != NULL; e_ptr = iter.get_next())
			e_ptr->finalize();
	}
}

void ScaffNode::clear_edges()
{
	outArcs.clr_list();
	inArcs.clr_list();
}

bool ScaffNode::isInternal()
{
	return (inArcs.get_size() <= 1 && outArcs.get_size() <= 1);
}

bool ScaffNode::isStart()
{
	return (outArcs.get_size() == 1);
}

bool ScaffNode::isEnd()
{
	return (inArcs.get_size() == 1);
}

bool ScaffNode::isEligible(ScaffNode * scaffolds)
{
	if(!this->isStart())
		return false;

	if(!this->isEnd())
		return true;	// otherwise the node has exactly one incoming edge

	LinkedIter<ScaffEdge> iter(inArcs);
	ScaffEdge * se = iter.get_first();
	return !(scaffolds[se->toContig].isStart());
}

bool ScaffNode::isLonger(int minLength)
{
	return (totalLen > minLength);
}
