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

#include <iostream>
#include <cassert>

#include "DnaRead.h"
#include "Naive.h"

using namespace std;

char rev_alphabet(char ch)
{
	switch(ch)
	{
		case 'A':
			return 'T';
		case 'T':
			return 'A';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
	}
	return 'N';
}

int DnaRead::PHRED = 64;

DnaRead::DnaRead()
{
	length = 0;
	seq = NULL;
	quals = NULL;
	rev = false;
	id = -1;
}

DnaRead::DnaRead(int max_length)
{
	length = max_length;
	seq = new char[max_length];
	quals = new char[max_length];
	seq[0] = '-';
	quals[0] = (char) (PHRED + 1);
	rev = false;
	id = -1;
}

DnaRead::~DnaRead()
{
	if(seq != NULL)
		delete [] seq;
	if(quals != NULL)
		delete [] quals;
}

void DnaRead::detach()
{
	seq = NULL;
	quals = NULL;
}

void DnaRead::replicate(const DnaRead & other)
{
	for(int i=0; i<other.length; i++)
		seq[i] = other.seq[i];

	length = other.length;
	id = other.id;
	rev= other.rev;
}

void DnaRead::replicate_with_quals(const DnaRead & other)
{
	for(int i=0; i<other.length; i++)
	{
		seq[i] = other.seq[i];
		quals[i] = other.quals[i];
	}

	length = other.length;
	id = other.id;
	rev= other.rev;
}

void DnaRead::replicate_partial(const DnaRead & other, int start, int end)
{
	int j = 1;
	for(int i=start; i<end; i++)
		seq[j++] = other.seq[i];

	length = j;
	id = other.id;
	rev= other.rev;
}

int DnaRead::copy_seq(const DnaRead & other, char * arr)
{
	seq = arr;
	for(int i=0; i<other.length; i++)
		seq[i] = other.seq[i];

	length = other.length;
	return length;
}

int DnaRead::copy_quals(const DnaRead & other, char * arr)
{
	assert(length == other.length);		// to make sure the sequence and quality values match
	quals = arr;

	for(int i=0; i<other.length; i++)
		quals[i] = other.quals[i];

	return length;
}

void DnaRead::revcomp()
{
	int len = length;
	int j = len-1;
	int mid = int(len/2) + 1;

	char tmp;
	char sc;
	for(int i=1; i<mid; i++)
	{
		tmp = seq[i];
		seq[i] = rev_alphabet( seq[j] );
		seq[j] = rev_alphabet( tmp );

		if(quals != NULL)
		{
			sc = quals[i];
			quals[i] = quals[j];
			quals[j] = sc;
		}

		j -= 1;
	}

	if(rev)
		rev = false;
	else
		rev = true;
}

void DnaRead::revcomp(bool state)
{
	if(rev != state)	// if not the right state reverse complement
	{
		int len = length;
		int j = len-1;
		int mid = int(len/2) + 1;

		char tmp;
		char sc;
		for(int i=1; i<mid; i++)
		{
			tmp = seq[i];
			seq[i] = rev_alphabet( seq[j] );
			seq[j] = rev_alphabet( tmp );

			if(quals != NULL)
			{
				sc = quals[i];
				quals[i] = quals[j];
				quals[j] = sc;
			}

			j -= 1;
		}

		if(rev)
			rev = false;
		else
			rev = true;
	}
}

bool DnaRead::trim(int trim1, int trim2, int droppoint)
{
	int len = length - trim1 - trim2;
	if(len < droppoint)
	{
		if(length < droppoint)	// never make a read larger than it is
			droppoint = length;

		for(int i=1; i<droppoint; i++)
		{
			seq[i] = 'N';
			quals[i] = (char) (PHRED + 1);
		}
		length = droppoint;
		return false;
	}

	int j = 1;
	for(int i=trim1+1; i<length-trim2; i++)
	{
		seq[j] = seq[i];
		quals[j++] = quals[i];
	}
	length = j;
	return true;
}

bool DnaRead::trim(int winsize, int droppoint, Naive * nb, double threshold)
{
	double total = 0.0;
	int trim1 = 0;
	int trim2 = 0;
	int lastbad = 0;

	for(int i=1; i<length; i++)
	{
		total += nb->log_correct[getq(i)];	// returns the negative log of the probability of the base being correct

		if(i > winsize)
			total -= nb->log_correct[getq(i - winsize)];

		if(total > threshold)
		{
			trim1 = i;
			lastbad = i;
		}
		if(i - lastbad >= droppoint)
			break;
	}

    lastbad = 0;
	total = 0.0;
	for(int i=1; i<length; i++)
	{
		total += nb->log_correct[getq(length - i)];

		if(i > winsize)
			total -= nb->log_correct[getq(length - i + winsize)];

		if(total > threshold)
		{
			trim2 = i;
			lastbad = i;
		}
		if(i - lastbad >= droppoint)
			break;
	}

	return trim(trim1, trim2, droppoint);
}

void DnaRead::write_fastq(ofstream & newReads, char * s, bool rc)
{
	revcomp(rc);	// make sure it's the right strand

	newReads << s << "\n";
	for(int j=1; j<length; j++)
		newReads << seq[j];
	newReads << "\n+\n";

	for(int j=1; j<length; j++)
		newReads << quals[j];

	newReads << "\n";
}

void DnaRead::write_fastq(ofstream & newReads, string & s, bool rc)
{
	revcomp(rc);	// make sure it's the right strand
	string s1 = s.substr(1);
	newReads << "@" << s1 << "\n";

	for(int j=1; j<length; j++)
		newReads << seq[j];
	newReads << "\n+\n";

	for(int j=1; j<length; j++)
		newReads << quals[j];

	newReads << "\n";
}

void DnaRead::write_fastq(string & s, bool rc)
{
	revcomp(rc);	// make sure it's the right strand
	string s1 = s.substr(1);
	cout << "@" << s1 << "\n";

	for(int j=1; j<length; j++)
		cout << seq[j];
	cout << "\n+\n";

	for(int j=1; j<length; j++)
		cout << quals[j];

	cout << "\n";
}

void DnaRead::write_seq(ofstream & outFile, int st, int end, int & bytes, int textwidth)
{
	for(int j=st; j<=end; j++)
	{
		outFile << seq[j];
		bytes++;
		if(bytes%textwidth == 0) outFile << "\n";
	}
}

bool DnaRead::get_kmer(int pos, int size, char * chars)
{
	pos = pos - (size/2);
	if(pos+size >= length)
		pos = length - size - 1;

	int i = 0;
	for(int j=pos; j<size+pos; j++)
	{
		char ch = seq[j];
		if(ch == 'N')
			return false;
		chars[i++] = ch;
	}
	return true;
}

void DnaRead::set_length(int len)
{
	length = len;
}
