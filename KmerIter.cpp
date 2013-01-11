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

#include "KmerNode.h"
#include "KmerList.h"
#include "KmerIter.h"

#ifndef NULL
	#define NULL 0
#endif

KmerIter::KmerIter(KmerList * L)
{
	if(L == NULL)
		throw "KmerIter: this list is null";
	list = L;
	cur = NULL;
}

KmerIter::~KmerIter()
{
	// nothing to do
}

KmerNode * KmerIter::get_first()
{
	cur = list->get_root();
	return cur;
}

KmerNode * KmerIter::get_cur()
{
	return cur;
}

KmerNode * KmerIter::get_next()
{
	if( cur != NULL )
		cur = cur->next;
	return cur;
}

void KmerIter::change_list(KmerList * L)
{
	if(L == NULL)
	{
		throw "KmerIter: List is null";
	}
	list = L;
	cur = NULL;
}

