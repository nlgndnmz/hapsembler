/*****************************************************************************
    $Author: Nilgun Donmez $
    $Date: 2011-05-05 18:58:27 -0400 (Thu, 05 May 2011) $

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

#ifndef NAIVE_H
#define NAIVE_H

class Naive
{
	public:

		Naive(int, int, int, double);
		Naive(const Naive &);
		~Naive();

		inline int ret_post(double w) { return int(w * num_steps); }
		void set_prior(double);

		double * log_correct;
		double * log_error;

		double ** log_same;
		double ** log_diff;

		double ** prior_same;
		double ** prior_diff;
		double ** not_prior;

		double priorB;
		double priorNotB;
		int probB;
		int num_steps;

		double * p_correct;
		double * p_error;
		double * prior;
		int max_qual;
};

#endif // NAIVE_H
