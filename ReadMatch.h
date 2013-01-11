/*****************************************************************************
    $Author: Nilgun Donmez $
    $Date: 2011-09-30 16:49:43 -0400 (Fri, 30 Sep 2011) $

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

#ifndef READ_MATCH_H
#define READ_MATCH_H

class ReadMatch
{
	public:

		ReadMatch(int);
		~ReadMatch();

		void reset();
		bool enter_hit(int);
		void compress(double, int, int, int, bool);
		bool get_next_pos(int &, double, int &, int &, int, int, int &, int);

		short int num_hits;

		short int * hits;
		short int * counts;

		int max_hits;
};

#endif // READ_MATCH_H
