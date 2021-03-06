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

#ifndef KMER_HASH_H
#define KMER_HASH_H

class KmerNode;
class KmerList;

class KmerHash
{
	public:
		KmerHash(int, long int, int);
		~KmerHash();

		KmerList * kmer_list;

		void add_entry(int, int, short int);
        bool search(int, int);

        int max_size;

	private:
		KmerNode * nodes;
		int hash_size;
		long int num_nodes;
		long int node_counter;

};

#endif // KMER_HASH_H

