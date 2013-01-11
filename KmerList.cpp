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

#ifndef NULL
	#define NULL 0
#endif

KmerList::KmerList()
{
	root = NULL;
	size = 0;
}

KmerList::~KmerList()
{
	// nothing to do
}

void KmerList::add_entry(int new_entry, short int position)
{
    KmerNode * newNode;
    newNode = new KmerNode();

    if (newNode == NULL)
    	throw "KmerList: Could not allocate a new node";

	newNode->data = new_entry;
	newNode->pos = position;

    if(root == NULL)    		// list was empty
        newNode->next = NULL;
    else
		newNode->next = root;
    root = newNode;				// the node becomes the new root

    size += 1;
}

void KmerList::clr_list()
{
	clr_list(root);
	size = 0;
}

void KmerList::clr_list(KmerNode * theNode)
{
	KmerNode * nextNode;
	while(theNode != NULL)
	{
		nextNode = theNode->next;
		delete theNode;
		theNode = nextNode;
	}
}

KmerNode * KmerList::get_root()
{
    return root;
}

int KmerList::get_size()
{
	return size;
}
