/*****************************************************************************

	Part of Hapsembler package. See the README file for more information.
    Copyright (C) 2011-2013,  Nilgun Donmez <nild@cs.toronto.edu>

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
#include "KmerHash.h"

#ifndef NULL
	#define NULL 0
#endif

KmerHash::KmerHash(int size, long int num_expected_nodes, int max_kmer_duplicity)
{
	hash_size = size;
	kmer_list = new KmerList[hash_size];

	num_nodes = num_expected_nodes + 5;		// so as not to complain for a small disagreement
	nodes = new KmerNode[num_nodes];
	node_counter = 0;
	max_size = max_kmer_duplicity;
}

KmerHash::~KmerHash()
{
	delete [] nodes;
	delete [] kmer_list;
}

void KmerHash::add_entry(int kmer, int id, short int position)
{
	if(node_counter >= num_nodes)
		throw "KmerHash: Maximum node capacity is reached!";

	if(kmer_list[kmer].size == max_size)	// reached the limit, add no more!
		return;

	KmerNode * newNode = &nodes[node_counter++];
	newNode->data = id;
	newNode->pos = position;

    if(kmer_list[kmer].root == NULL)    // list was empty
    	newNode->next = NULL;
    else
        newNode->next = kmer_list[kmer].root;

    kmer_list[kmer].root = newNode;
    kmer_list[kmer].size += 1;
}

bool KmerHash::search(int kmer, int id)
{
	KmerNode * nd = kmer_list[kmer].root;
	while(nd != NULL)
	{
		if(nd->data == id)
			return true;
		nd = nd->next;
	}
	return false;
}
