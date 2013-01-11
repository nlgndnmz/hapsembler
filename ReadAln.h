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

#ifndef READ_ALN_H
#define READ_ALN_H

class ReadAln
{
	public:

		ReadAln();
		~ReadAln();

		void set(int);

		char * read1;
		char * read0;
		short int * qual1;
		short int * qual0;
		int pos;
		int start;
		int weight;
		int longest_gap;
		int diff;
};

#endif // READ_ALN_H
