/*****************************************************************************
    $Author: Nilgun Donmez $
    $Date: 2012-10-05 17:38:02 -0400 (Fri, 05 Oct 2012) $

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

#ifndef DNA_READ_H
#define DNA_READ_H

#include <fstream>
#include <string>

class Naive;

class DnaRead
{
	public:
		DnaRead();
		explicit DnaRead(int);
		~DnaRead();

		void detach();

		void revcomp();
		void revcomp(bool);

		static void set_phred(int phred_offset) { PHRED = phred_offset; };

		void write_fastq(std::ofstream &, char *, bool = false);
		void write_fastq(std::ofstream &, std::string &, bool = false);
		void write_fastq(std::string &, bool = false);
		void write_seq(std::ofstream &, int, int, int &, int);

		void replicate(const DnaRead &);
		void replicate_with_quals(const DnaRead &);
		void replicate_partial(const DnaRead &, int, int);

		int copy_quals(const DnaRead &, char *);
		int copy_seq(const DnaRead &, char *);

		inline int getq(int i) { return ((int)(quals[i]) - PHRED); };
		bool trim(int, int, Naive *, double);
		bool trim(int, int, int);

		void set_length(int);

		bool get_kmer(int, int, char *);

		char * seq;
		char * quals;

		int length;
		bool rev;
		int id;

		static int PHRED;
};

#endif // DNA_READ_H

