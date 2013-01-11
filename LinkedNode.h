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

#ifndef LINKED_NODE_H
#define LINKED_NODE_H

#include <iostream>

template <class TYPE>
class LinkedNode
{
    public:
		void * operator new(size_t);
		void operator delete(void *);

        explicit LinkedNode(TYPE *);
        LinkedNode(TYPE *, char);
        LinkedNode();
        ~LinkedNode();

        TYPE * data;
        LinkedNode<TYPE> * next;
        char flg;	// flag

		static void shrink();
		static void purge();

	private:

		bool alive;

		static void init();

		static LinkedNode ** pool;

		static int poolnum;
		static int chunkctr;

		const static int chunksize = 1048576;	// 1024x1024
		const static int poolsize = 1048576;
};

// define the static variables
template <class TYPE> LinkedNode<TYPE> ** LinkedNode<TYPE>::pool = 0;
template <class TYPE> int LinkedNode<TYPE>::poolnum = 0;
template <class TYPE> int LinkedNode<TYPE>::chunkctr = 0;

template <class TYPE>
void * LinkedNode<TYPE>::operator new(size_t sz)
{
	if(pool == 0)
		init();

	if(chunkctr < chunksize)	// i.e. if there is still some space in the current chunk
	{
		chunkctr += 1;
		return &((pool[poolnum])[chunkctr-1]);
	}
	else if((poolnum+1) < poolsize)		// then allocate the next chunk
	{
		pool[++poolnum] = new LinkedNode<TYPE>[chunksize];
		chunkctr = 1;
		return &((pool[poolnum])[chunkctr-1]);
	}
	else	// out of memory, throw an exception
		throw "LinkedNode pool is out of memory!";
}

template <class TYPE>
void LinkedNode<TYPE>::operator delete(void * p)
{
	// don't do anything
}

template <class TYPE>
void LinkedNode<TYPE>::shrink()
{
	//std::cout << "Shrinking LinkedNode pool" << std::endl;
	if(pool != 0)
	{
		for(int i=0; i<poolnum; i++)
		{
			bool alldead = true;
			for(int j=0; j<chunksize; j++)
			{
				if(pool[i][j].alive)
					alldead = false;
			}
			if(alldead)
			{
				delete [] pool[i];
				pool[i] = 0;
			}
		}
	}
}

template <class TYPE>
void LinkedNode<TYPE>::purge()
{
	//std::cout << "LinkedNode pool is being purged" << std::endl;
	if(pool != 0)
	{
		for(int i=0; i<=poolnum; i++)
			delete [] pool[i];

		delete [] pool;
	}
	pool = 0;
	poolnum = 0;
	chunkctr = 0;
}

template <class TYPE>
void LinkedNode<TYPE>::init()
{
	//std::cout << "LinkedNode pool is being created" << std::endl;
	pool = new LinkedNode<TYPE>*[poolsize];
	pool[0] = new LinkedNode<TYPE>[chunksize];
	chunkctr = 0;
}

template<class TYPE>
LinkedNode<TYPE>::LinkedNode()
{
   data = 0;
   flg = 0;
   alive = true;
}

template<class TYPE>
LinkedNode<TYPE>::LinkedNode(TYPE * the_data)
{
   data = the_data;
   flg = 0;
   alive = true;
}

template<class TYPE>
LinkedNode<TYPE>::LinkedNode(TYPE * the_data, char w)
{
   data = the_data;
   flg = w;
   alive = true;
}

template<class TYPE>
LinkedNode<TYPE>::~LinkedNode()
{
    alive = false;
}

#endif
