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

#ifndef LINKED_ITER_H
#define LINKED_ITER_H

#include "LinkedNode.h"
#include "LinkedList.h"

template <class TYPE>
class LinkedIter
{
public:
	explicit LinkedIter(LinkedList<TYPE> &);
	LinkedIter();
	~LinkedIter();

	TYPE * get_first();
	TYPE * get_next();
	TYPE * del_cur();
	TYPE * del_cur(bool);
	TYPE * get_cur();
	TYPE * insert_new(TYPE *);

	int get_flg();
	void change_list(LinkedList<TYPE> &);
	int get_size();

private:
	LinkedNode<TYPE> * prev;
	LinkedNode<TYPE> * cur;
	LinkedList<TYPE> * list;

};

template<class TYPE>
LinkedIter<TYPE>::LinkedIter(LinkedList<TYPE> & L)
{
	list = &L;
	cur = NULL;
	prev = NULL;
}

template<class TYPE>
LinkedIter<TYPE>::LinkedIter()
{
	list = NULL;
	cur = NULL;
	prev = NULL;
}

template<class TYPE>
LinkedIter<TYPE>::~LinkedIter()
{
	// nothing to do
}

template<class TYPE>
int LinkedIter<TYPE>::get_flg()
{
	assert(cur != NULL);
	return (int) (cur->flg);
}

template<class TYPE>
int LinkedIter<TYPE>::get_size()
{
	assert(list != NULL);
	return list->get_size();
}

template<class TYPE>
TYPE * LinkedIter<TYPE>::get_first()
{
	cur = list->get_root();
	prev = NULL;
	if(cur != NULL)
		return cur->data;

	return NULL;
}

template<class TYPE>
TYPE * LinkedIter<TYPE>::get_cur()
{
	if( cur != NULL )
		return cur->data;

	return NULL;
}

template<class TYPE>
TYPE * LinkedIter<TYPE>::get_next()
{
	prev = cur;
	if( cur != NULL )
	{
		cur = cur->next;
		if(cur != NULL)
			return cur->data;
	}
	return NULL;
}

template<class TYPE>
void LinkedIter<TYPE>::change_list(LinkedList<TYPE> & L)
{
	list = &L;
	cur = NULL;
	prev = NULL;
}

template<class TYPE>
TYPE * LinkedIter<TYPE>::del_cur(bool override)
{
	assert(cur != NULL);

	LinkedNode<TYPE> * temp = cur->next;
	list->del_entry(prev, override);	// prev remains as before
	cur = temp;
	if(cur != NULL)
		return cur->data;

	return NULL;
}

template<class TYPE>
TYPE * LinkedIter<TYPE>::del_cur()
{
	assert(cur != NULL);

	LinkedNode<TYPE> * temp = cur->next;
	list->del_entry(prev);
	cur = temp;
	if(cur != NULL)
		return cur->data;

	return NULL;
}

template<class TYPE>
TYPE * LinkedIter<TYPE>::insert_new(TYPE * the_data)
{
	LinkedNode<TYPE> * temp = cur->next;
	prev = list->insert_entry(the_data, prev);
	cur = temp;
	if(cur != NULL)
		return cur->data;

	return NULL;
}

#endif
