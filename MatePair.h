/*****************************************************************************
    $Author: Nilgun Donmez $
    $Date: 2012-01-26 20:21:34 -0500 (Thu, 26 Jan 2012) $

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

#ifndef MATE_PAIR_H
#define MATE_PAIR_H

#include <fstream>

#include "GraphDef.h"
#include "BiNode.h"
#include "MateEdge.h"

class MatePair : public BiNode<MateEdge>
{
	public:

		MatePair();
		~MatePair();

		void pickle(std::ofstream &);
        void unpickle(std::ifstream &, int, int, MatePair *);

        bool isEligible();

		short int lib;

        const static ClassID data_id = MATE_PAIR;
        static double FLEX;
};

double MatePair::FLEX = 1.0001;

MatePair::MatePair()
{
	lib = 0;
}

MatePair::~MatePair()
{
	// base class handles the rest
}

// pickle the node to the file given
void MatePair::pickle(std::ofstream & fh)
{
	int num = get_degree();
	fh << (int)data_id <<" "<< num <<" "<< id <<" "<< len <<" "<< flag <<" "<< lib <<" "<< FLEX <<"\n";
	BiNode<MateEdge>::pickle(fh);
}

// unpickle the node from file
void MatePair::unpickle(std::ifstream & fh, int num, int id, MatePair * pairs)
{
	this->id = id;
	fh >> len >> flag >> lib >> FLEX;

	for(int i=0; i<num; i++)
	{
        MateEdge * me = new MateEdge();
		me->unpickle(fh);
		int k = (me->toRead[0] == id) ? 1 : 0;		// nodes should not have self-loops
		if(id < me->toRead[k] || (id == me->toRead[k] && k==1))
		{
            add_arc(me, k);
            pairs[me->toRead[k]].add_arc(me, rev(k));
		}
		else
            delete me;
	}
}

bool MatePair::isEligible()
{
	return (lib != 0);
}

#endif // MATE_PAIR_H
