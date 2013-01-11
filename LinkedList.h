/*****************************************************************************
    $Author: Nilgun Donmez $
    $Date: 2011-06-08 13:43:52 -0400 (Wed, 08 Jun 2011) $

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

#ifndef LINKED_LIST_H
#define LINKED_LIST_H

#include <cassert>
#include "LinkedNode.h"

#ifndef NULL
	#define NULL 0
#endif

template <class TYPE>
class LinkedList
{
    public:
        explicit LinkedList(bool);
        LinkedList();
        ~LinkedList();

        void add_entry(TYPE *);
        void add_entry(TYPE *, char);
        void add_sorted(TYPE *, char);

        void sort();
        void sort2();
        void check();

        void insert_entry(TYPE *);
        void insert_entry(TYPE *, char);
        LinkedNode<TYPE>* insert_entry(TYPE *, LinkedNode<TYPE> *);
        void del_entry(LinkedNode<TYPE> *);
        void del_entry(LinkedNode<TYPE> *, bool);

        void append_list(LinkedList<TYPE> *);	// sort of deep copy but appends to the tail

        void set_delete();
        void set_no_delete();

        LinkedNode<TYPE> * get_root();
        LinkedNode<TYPE> * get_tail();

        int get_size();
        void clr_list();

    private:

		void clr_list(LinkedNode<TYPE> *);
        void heapify(LinkedNode<TYPE> **, int, int);
        void heapify2(LinkedNode<TYPE> **, int, int);

        LinkedNode<TYPE> * root;
        LinkedNode<TYPE> * tail;

        bool del_data;    // decides if the actual data should be deleted or not
        int num_loaded;
};

template<class TYPE>
LinkedList<TYPE>::LinkedList(bool del_flag)
{
	root = NULL;
	tail = NULL;
	num_loaded = 0;
	del_data = del_flag;    // force the client to decide that in advance
}

template<class TYPE>
LinkedList<TYPE>::LinkedList()
{
	root = NULL;
	tail = NULL;
	num_loaded = 0;
	del_data = true;		// default behaviour but user may override later
}

template<class TYPE>
LinkedList<TYPE>::~LinkedList()
{
    clr_list(root);
	root = NULL;
	tail = NULL;
	num_loaded = 0;
}

template<class TYPE>
void LinkedList<TYPE>::set_delete()
{
	del_data = true;
}

template<class TYPE>
void LinkedList<TYPE>::set_no_delete()
{
	del_data = false;
}

// clr_list - The actual TYPE data is not deleted if del_data is 0
template<class TYPE>
void LinkedList<TYPE>::clr_list(LinkedNode<TYPE> * theNode)
{
	LinkedNode<TYPE> * nextNode;
	while(theNode != NULL)
	{
		nextNode = theNode->next;
		if(del_data)    // delete data if requested
		{
		   if(theNode->data != NULL)    // check if NULL
		   {
		       delete theNode->data;
		       theNode->data = NULL;
		   }
		}
		delete theNode;
		theNode = nextNode;
		num_loaded--;
	}
}

template<class TYPE>
void LinkedList<TYPE>::clr_list()
{
	LinkedNode<TYPE> * nextNode;
	LinkedNode<TYPE> * theNode = root;
	while(theNode != NULL)
	{
		nextNode = theNode->next;
		if(del_data)    // delete data if requested
		{
		   if(theNode->data != NULL)    // check if NULL
		   {
		       delete theNode->data;
		       theNode->data = NULL;
		   }
		}
		delete theNode;
		theNode = nextNode;
	}
	root = NULL;
	tail = NULL;
	num_loaded = 0;
}

// Adds the entry at the beginning of the list
template<class TYPE>
void LinkedList<TYPE>::insert_entry(TYPE *new_entry)
{
    assert(new_entry != NULL);

    LinkedNode<TYPE> *newNode;
    newNode = new LinkedNode<TYPE>(new_entry);

    if(newNode == NULL)
    {	throw "LinkedList: Could not allocate a new node";	}

    if(root == NULL)    // list was empty
    {
        root = newNode;
        tail = newNode;
        newNode->next = NULL;
    }
    else
    {
        newNode->next = root;
        root = newNode;		// new root
    }

    num_loaded++;
}

// Adds the entry at the beginning of the list
template<class TYPE>
void LinkedList<TYPE>::insert_entry(TYPE *new_entry, char which)
{
    assert(new_entry != NULL);

    LinkedNode<TYPE> *newNode;
    newNode = new LinkedNode<TYPE>(new_entry, which);

    if(newNode == NULL)
    {	throw "LinkedList: Could not allocate a new node";	}

    if(root == NULL)    // list was empty
    {
        root = newNode;
        tail = newNode;
        newNode->next = NULL;
    }
    else
    {
        newNode->next = root;
        root = newNode;		// new root
    }

    num_loaded++;
}


// Adds the entry at the after the given node
template<class TYPE>
LinkedNode<TYPE>* LinkedList<TYPE>::insert_entry(TYPE *new_entry, LinkedNode<TYPE> * nd)
{
    assert(new_entry != NULL);
    assert(nd != NULL);

    LinkedNode<TYPE> *newNode;
    newNode = new LinkedNode<TYPE>(new_entry);

    if(newNode == NULL)
    {	throw "LinkedList: Could not allocate a new node";	}

    LinkedNode<TYPE> * nextNode = nd->next;

    if(nextNode == NULL)    // nd was the tail
    {
        tail = newNode;
        newNode->next = NULL;
        nd->next = newNode;
    }
    else
    {
       newNode->next = nextNode;
       nd->next = newNode;
    }

    num_loaded++;
    return newNode;
}

// add_entry - Adds at the end of the list
template<class TYPE>
void LinkedList<TYPE>::add_entry(TYPE *new_entry)
{
    assert(new_entry != NULL);

    LinkedNode<TYPE> *newNode;
    newNode = new LinkedNode<TYPE>(new_entry);

    if(newNode == NULL)
    {	throw "LinkedList: Could not allocate a new node";	}

    newNode->next = NULL;
    if(root == NULL)    // list was empty
    {
        root = newNode;
        tail = newNode;
    }
    else
    {
        tail->next = newNode;
        tail = newNode;		// new tail
    }

    num_loaded++;
}

// add_entry - Adds at the end of the list
template<class TYPE>
void LinkedList<TYPE>::add_entry(TYPE *new_entry, char which)
{
    assert(new_entry != NULL);

    LinkedNode<TYPE> *newNode;
    newNode = new LinkedNode<TYPE>(new_entry, which);

    if(newNode == NULL)
    {	throw "LinkedList: Could not allocate a new node";	}

    newNode->next = NULL;
    if(root == NULL)    // list was empty
    {
        root = newNode;
        tail = newNode;
    }
    else
    {
        tail->next = newNode;
        tail = newNode;		// new tail
    }

    num_loaded++;
}

template<class TYPE>
void LinkedList<TYPE>::append_list(LinkedList<TYPE> * other_list)
{
	LinkedNode<TYPE> * nd = other_list->root;
	while(nd != NULL)
	{
		assert(nd->data != NULL);

		TYPE * t = new TYPE(*(nd->data));
		this->add_entry(t);
		nd = nd->next;
	}
}

template<class TYPE>
void LinkedList<TYPE>::check()
{
	LinkedNode<TYPE> * nd = root;
	int counter = 0;
	while(nd != NULL)
	{
		counter++;
		assert(counter <= num_loaded);
		nd = nd->next;
	}
}

// Sorts the list
template<class TYPE>
void LinkedList<TYPE>::sort()
{
    if(root == NULL)    // if the list is empty nothing to sort
        return;

	LinkedNode<TYPE> ** maxheap = new LinkedNode<TYPE>*[num_loaded];

	// populate the heap
	LinkedNode<TYPE> * nd = root;
	int index = 0;
	while(nd != NULL)
	{
		maxheap[index++] = nd;
		nd = nd->next;
	}

	assert(index == num_loaded);

	int start = int((index - 2) / 2.0);
	while(start >= 0)
	{
		heapify(maxheap, start, index);
		start--;
	}

	index -= 1;
	while(index > 0)
	{
		LinkedNode<TYPE> * tmp = maxheap[index];
		maxheap[index] = maxheap[0];
		maxheap[0] = tmp;
		heapify(maxheap, 0, index--);
	}

	root = maxheap[0];
	for(int i=0; i<num_loaded-1; i++)
		maxheap[i]->next = maxheap[i+1];

	tail = maxheap[num_loaded-1];
	tail->next = NULL;

	delete [] maxheap;
}

template<class TYPE>
void LinkedList<TYPE>::heapify(LinkedNode<TYPE> ** maxheap, int head, int end)
{
	LinkedNode<TYPE> * head_ptr = maxheap[head];

	while(head * 2 + 1 < end)
	{
		int left_child = head * 2 + 1;
		int right_child = left_child + 1;

		LinkedNode<TYPE> * left = maxheap[left_child];
		LinkedNode<TYPE> * child_ptr = left;
		int child = left_child;

		if(right_child < end)
		{
			LinkedNode<TYPE> * right = maxheap[right_child];
			if(left->data->len[(int)left->flg] < right->data->len[(int)right->flg])
			{
				child_ptr = right;
				child = right_child;
			}
		}

		if( head_ptr->data->len[(int)head_ptr->flg] < child_ptr->data->len[(int)child_ptr->flg] )
		{
			maxheap[head] = maxheap[child];
			maxheap[child] = head_ptr;
			head = child;
		}
		else
			return;
	}
}

// Sorts the list
template<class TYPE>
void LinkedList<TYPE>::sort2()
{
    if(root == NULL)    // if the list is empty nothing to sort
        return;

	LinkedNode<TYPE> ** maxheap = new LinkedNode<TYPE>*[num_loaded];

	// populate the heap
	LinkedNode<TYPE> * nd = root;
	int index = 0;
	while(nd != NULL)
	{
		maxheap[index++] = nd;
		nd = nd->next;
	}

	assert(index == num_loaded);

	int start = int((index - 2) / 2.0);
	while(start >= 0)
	{
		heapify(maxheap, start, index);
		start--;
	}

	index -= 1;
	while(index > 0)
	{
		LinkedNode<TYPE> * tmp = maxheap[index];
		maxheap[index] = maxheap[0];
		maxheap[0] = tmp;
		heapify(maxheap, 0, index--);
	}

	root = maxheap[0];
	for(int i=0; i<num_loaded-1; i++)
		maxheap[i]->next = maxheap[i+1];

	tail = maxheap[num_loaded-1];
	tail->next = NULL;

	delete [] maxheap;
}

template<class TYPE>
void LinkedList<TYPE>::heapify2(LinkedNode<TYPE> ** maxheap, int head, int end)
{
	LinkedNode<TYPE> * head_ptr = maxheap[head];

	while(head * 2 + 1 < end)
	{
		int left_child = head * 2 + 1;
		int right_child = left_child + 1;

		LinkedNode<TYPE> * left = maxheap[left_child];
		LinkedNode<TYPE> * child_ptr = left;
		int child = left_child;

		if(right_child < end)
		{
			LinkedNode<TYPE> * right = maxheap[right_child];
			if(left->data->len[0] < right->data->len[0])
			{
				child_ptr = right;
				child = right_child;
			}
		}

		if( head_ptr->data->len[0] < child_ptr->data->len[0] )
		{
			maxheap[head] = maxheap[child];
			maxheap[child] = head_ptr;
			head = child;
		}
		else
			return;
	}
}

// add_sorted - Adds the entry sorted by the key
template<class TYPE>
void LinkedList<TYPE>::add_sorted(TYPE *new_entry, char which)
{
    assert(new_entry != NULL);

    LinkedNode<TYPE> * newNode;
    newNode = new LinkedNode<TYPE>(new_entry, which);

    short int k = new_entry->len[(int)which];

    if(newNode == NULL)
    {	throw "LinkedList: Could not allocate a new node";	}

    if(root == NULL)    // if the list is empty
    {
        root = newNode;
        tail = newNode;
        newNode->next = NULL;
    }
    else
    {
	    LinkedNode<TYPE> * nd = root;
	    LinkedNode<TYPE> * nd_prev = NULL;
	    while(nd != NULL && nd->data->len[(int) nd->flg] < k)
	    {
	    	nd_prev = nd;
		    nd = nd->next;
	    }

	    // either nd is NULL or nd's key is bigger
	    // so we need to insert the new node before nd
    	newNode->next = nd;

	    if(nd != NULL)
	    {
	    	if(nd_prev != NULL)
				nd_prev->next = newNode;
			else
				root = newNode;
		}
		else	// otherwise, it means we are adding to the tail
		{
			tail->next = newNode;
			tail = newNode;		// new tail
		}
	}

    num_loaded++;
}

// Deletes the entry next to the one pointed by cur, overrides the del_data field by doing real deletion
template<class TYPE>
void LinkedList<TYPE>::del_entry(LinkedNode<TYPE> * prevNode, bool override)
{
	LinkedNode<TYPE> * cur;
	if(prevNode == NULL)
		cur = root;
	else
		cur = prevNode->next;

	assert(cur != NULL);

	LinkedNode<TYPE> * nextNode = cur->next;

	// unlink the current node
	if(prevNode != NULL)
		prevNode->next = nextNode;
	else    // if the current node was root
		root = nextNode;

	if(nextNode == NULL)
		tail = prevNode;

	if(override && cur->data != NULL)
	{
		delete cur->data;
		cur->data = NULL;
	}

	delete cur;
	cur = NULL;

	num_loaded--;
}

// del_cur_entry
template<class TYPE>
void LinkedList<TYPE>::del_entry(LinkedNode<TYPE> * prevNode)
{
	LinkedNode<TYPE> * cur;
	if(prevNode == NULL)
		cur = root;
	else
		cur = prevNode->next;

	assert(cur != NULL);

	LinkedNode<TYPE> * nextNode = cur->next;

	// unlink the current node
	if(prevNode != NULL)
		prevNode->next = nextNode;
	else    // if the current node was root
		root = nextNode;

	if(nextNode == NULL)
		tail = prevNode;

    if(del_data && cur->data != NULL)
    {
		delete cur->data;
		cur->data = NULL;
    }

    delete cur;
    cur = NULL;

    num_loaded--;
}

// get_first - gets the root
template<class TYPE>
LinkedNode<TYPE> * LinkedList<TYPE>::get_root()
{
    if(root != NULL)
    	return root;

    return NULL;
}

// get_last - gets the tail
template<class TYPE>
LinkedNode<TYPE> * LinkedList<TYPE>::get_tail()
{
    if(tail != NULL)
    	return tail;

    return NULL;
}

// get_num_loaded - gets the size of the list
template<class TYPE>
int LinkedList<TYPE>::get_size()
{
	return num_loaded;
}

#endif // LINKED_LIST_H

