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

#include <cassert>
#include <cmath>

#include "SmithWaterman.h"
#include "HapUtils.h"
#include "DnaRead.h"
#include "ReadAln.h"
#include "Naive.h"

using namespace std;

SmithWaterman::SmithWaterman(int max_read_size, int min_overlap, int match, int mismatch, int gap, double epsilon, int gap_score,
	int max_quality, int num_nucleotides, int max_steps, double prior)
{
	this->max_read_size = max_read_size;
	this->min_overlap = min_overlap;
	this->epsilon = epsilon;

	GAP = gap;
	MISMATCH = mismatch;
	MATCH = match;
	FORBIDDEN = GAP * max_read_size;

	matrix = 0;	// NULL

	this->gap_score = gap_score;
	this->NB = new Naive(max_quality, num_nucleotides, max_steps, prior);

	numCached = 0;
	maxCache = 0;
	theCache = 0;
}

SmithWaterman::SmithWaterman(int max_read_size, int min_overlap, int match, int mismatch, int gap, double epsilon)
{
	this->max_read_size = max_read_size;
	this->min_overlap = min_overlap;
	this->epsilon = epsilon;

	GAP = gap;
	MISMATCH = mismatch;
	MATCH = match;
	FORBIDDEN = GAP * max_read_size;

	this->NB = 0; // NULL
	this->set_matrix(max_read_size);
	this->set_cache(2);

	numCached = 0;
	maxCache = 2;
}

SmithWaterman::SmithWaterman(const SmithWaterman & sw)
{
	this->max_read_size = sw.max_read_size;
	this->min_overlap = sw.min_overlap;
	this->epsilon = sw.epsilon;

	this->GAP = sw.GAP;
	this->MISMATCH = sw.MISMATCH;
	this->MATCH = sw.MATCH;

	this->set_matrix(max_read_size);

	this->gap_score = sw.gap_score;
	this->NB = new Naive(*(sw.NB));

	numCached = 0;
	maxCache = 0;
}

SmithWaterman::~SmithWaterman()
{
	if(matrix != 0)
	{
		for(int i=0; i<max_read_size; i++)
			delete [] matrix[i];

		delete [] matrix;
	}
	if(NB != 0)
		delete NB;

	if(theCache != 0)
		delete [] theCache;
}

void SmithWaterman::set_matrix(int max_rs)
{
	max_read_size = max_rs;
	matrix = new int*[max_read_size];
	FORBIDDEN = GAP * max_read_size;

	for(int i=0; i<max_read_size; i++)
	{
		matrix[i] = new int[max_read_size];
		matrix[i][0] = 0;
	}

	for(int j=0; j<max_read_size; j++)
		matrix[0][j] = 0;
}

void SmithWaterman::set_cache(int size)
{
	theCache = new OvlCache[size];
	maxCache = size;
	numCached = 0;
}

void SmithWaterman::reset_cache()
{
	numCached = 0;
}

int SmithWaterman::get_numCached()
{
	return numCached;
}

OvlCache * SmithWaterman::get_cache()
{
	return theCache;
}

void SmithWaterman::flush_cache(ofstream & ovlFile)
{
	for(int i=0; i<numCached; i++)
	{
		ovlFile << theCache[i].IDs[0] << " " << theCache[i].IDs[1];
		if(theCache[i].strand)
			ovlFile << " 1 ";
		else
			ovlFile << " 0 ";

		for(int j=0; j<4; j++)
			ovlFile << theCache[i].pos[j] << " ";

		ovlFile << theCache[i].weight <<" "<< theCache[i].indels <<"\n";
	}
	numCached = 0;
}

int SmithWaterman::align(DnaRead & rd1, DnaRead & rd2, int k1, int k2, int width, int & pos, ReadAln & rdaln,
	bool swap_reads, bool isConsensr)
{
	char * str1 = rd1.seq;
	char * str2 = rd2.seq;

	int len1 = rd1.length;
	int len2 = rd2.length;

	int max_score = 0;
	int row, col=-1, max_row=0, max_col=0;

	int rowStart = k1;
	int colStart = k2;
	int colEnd;

	for(int i=0; i<len1; i++)
	{
		if(isConsensr)
			matrix[i][0] = FORBIDDEN;
		else
			matrix[i][0] = 0;

		matrix[i][len2-1] = FORBIDDEN;
	}

	for(int j=0; j<len2; j++)
	{
		matrix[0][j] = 0;
		matrix[len1-1][j] = FORBIDDEN;
	}

	if(k1 > 1)
		matrix[k1-1][1] = FORBIDDEN;

	for(row=rowStart; row<len1; row++)
	{
		colStart++;
		if(colStart >= len2)			// no need to go further, already processed the entire band
			break;

		colEnd = colStart + width;
		colEnd = (colEnd > len2) ? len2 : colEnd;
		col = (colStart > 0) ? colStart : 1;

		if(col > 1)
			matrix[row][col-1] = FORBIDDEN;

		for(; col<colEnd; col++)
		{
			int up, left, diagonal, d, better;

			d = (str1[row] == str2[col]) ? MATCH : MISMATCH;

			diagonal = matrix[row-1][col-1] + d;
			up = matrix[row-1][col];
			left = matrix[row][col-1];

			better = (up > left) ? up + GAP : left + GAP;
			matrix[row][col] = (better > diagonal) ? better : diagonal;

		}
		if(col < len2)
			matrix[row][col] = FORBIDDEN;	// we set the boundary to something very negative, so that this path will never be taken
	}

	assert(col > 0);

	for(int i=1; i<len1; i++)
	{
		if(matrix[i][len2-1] >= max_score)
		{
			max_score = matrix[i][len2-1];
			max_row = i;
			max_col = len2-1;
		}
	}

	for(int j=1; j<len2; j++)
	{
		if(matrix[len1-1][j] >= max_score)
		{
			max_score = matrix[len1-1][j];
			max_row = len1-1;
			max_col = j;
		}
	}

	col = max_col;
	row = max_row;

	char * line1;	// these are just namesakes
	char * line2;
	short int * score1;
	short int * score2;

	if(swap_reads)
	{
		line1 = rdaln.read1;
		line2 = rdaln.read0;

		score1 = rdaln.qual1;
		score2 = rdaln.qual0;
	}
	else
	{
		line1 = rdaln.read0;
		line2 = rdaln.read1;

		score1 = rdaln.qual0;
		score2 = rdaln.qual1;
	}

	int index = 0;

	double pB = 0.0;
	double pnotB = 0.0;

	int weight = 0;
	int insertions = 0;
	int deletions = 0;
	short int sc1 = 0;
	short int sc2 = 0;

	int longest_gap = 0;
	int running_gap = 0;
	char last_char = ' ';

	while(row > 0 && col > 0)
	{
		int d = (str1[row] == str2[col]) ? MATCH : MISMATCH;
		int score = matrix[row][col];
		char arr = 'L';

		if(score == (matrix[row-1][col-1] + d))
			arr = 'D';
		else if(score == (matrix[row-1][col] + GAP))
			arr = 'U';

		if(last_char == arr)
			running_gap++;
		else
		{
			if(running_gap > longest_gap)
				longest_gap = running_gap;
			running_gap = 0;
			if(arr != 'D')
				last_char = arr;
			else
				last_char = ' ';
		}

		switch(arr)
		{
			case 'D':

				sc1 = ((short int) rd1.quals[row]) - DnaRead::PHRED;
				sc2 = ((short int) rd2.quals[col]) - DnaRead::PHRED;

				if(str1[row] != str2[col] && str1[row] != 'N' && str2[col] != 'N')
				{
                    pB += NB->prior_diff[ sc1 ][ sc2 ];
                    pnotB += NB->not_prior[ sc1 ][ sc2 ];
                    weight++;
				}

				score1[index] = sc1;
				score2[index] = sc2;

				line1[index] = str1[row];
				line2[index++] = str2[col];

				row--; col--;
				break;

			case 'L':

				sc2 = ((short int) rd2.quals[col]) - DnaRead::PHRED;

				pB += NB->prior_diff[ gap_score ][ sc2 ];
				pnotB += NB->not_prior[ gap_score ][ sc2 ];

				score1[index] = gap_score;
				score2[index] = sc2;

				line1[index] = '-';
				line2[index++] = str2[col];

				weight++;
				deletions++;
				col--;
				break;

			case 'U':

				sc1 = ((short int) rd1.quals[row]) - DnaRead::PHRED;

				pB += NB->prior_diff[ sc1 ][ gap_score ];
				pnotB += NB->not_prior[ sc1 ][ gap_score ];

				score1[index] = sc1;
				score2[index] = gap_score;

				line1[index] = str1[row];
				line2[index++] = '-';

				weight++;
				insertions++;
				row--;
				break;
		}

	}

	int ovl_size = (max_row - row) + (max_col - col);
	ovl_size = int(ovl_size/2);

	sc1 = 0;
	sc2 = 0;

	if(weight > 0)
    {
        pB = pB / weight;
        pnotB = pnotB / weight;
    }

	pB += NB->priorB;
	pnotB += NB->priorNotB;

	rdaln.weight = NB->ret_post(1.0 / ( 1.0 + exp(pnotB - pB)));
	rdaln.longest_gap = longest_gap;
	rdaln.diff = weight;

	pos = col - row;	// at least one of them should be zero anyways, so this should give the correct offset in the case of consensr
	return (index-1);	// returns the last index of the overlap
}

// Special alignment routine for overlappr only
bool SmithWaterman::align(DnaRead & rd1, DnaRead & rd2, int k1, int k2, int width, bool strand, bool & isContained,
	bool nogap)
{
	int max_score = 0;
	int row, col=-1, max_row=0, max_col=0;

	int rowStart = k1;
	int colStart = k2;
	int colEnd;

	for(int i=0; i<rd1.length; i++)
	{
		matrix[i][0] = 0;
		matrix[i][rd2.length-1] = FORBIDDEN;
	}

	for(int j=0; j<rd2.length; j++)
	{
		matrix[0][j] = 0;
		matrix[rd1.length-1][j] = FORBIDDEN;
	}

	if(k1 > 1)
		matrix[k1-1][1] = FORBIDDEN;

	for(row=rowStart; row<rd1.length; row++)
	{
		colStart++;
		if(colStart >= rd2.length)			// no need to go further, already processed the entire band
			break;

		colEnd = colStart + width;
		colEnd = (colEnd > rd2.length) ? rd2.length : colEnd;
		col = (colStart > 0) ? colStart : 1;

		if(col > 1)
			matrix[row][col-1] = FORBIDDEN;

		for(; col<colEnd; col++)
		{
			int up, left, diagonal, d, better;

			d = (rd1.seq[row] == rd2.seq[col]) ? MATCH : MISMATCH;

			diagonal = matrix[row-1][col-1] + d;
			up = matrix[row-1][col];
			left = matrix[row][col-1];

			better = (up > left) ? up + GAP : left + GAP;
			matrix[row][col] = (better > diagonal) ? better : diagonal;

		}
		if(col < rd2.length)
			matrix[row][col] = FORBIDDEN;	// we set the boundary to something very negative, so that this path will never be taken
	}

	assert(col > 0);

	for(int i=1; i<rd1.length; i++)
	{
		if(matrix[i][rd2.length-1] >= max_score)
		{
			max_score = matrix[i][rd2.length-1];
			max_row = i;
			max_col = rd2.length-1;
		}
	}

	for(int j=1; j<rd2.length; j++)
	{
		if(matrix[rd1.length-1][j] >= max_score)
		{
			max_score = matrix[rd1.length-1][j];
			max_row = rd1.length-1;
			max_col = j;
		}
	}

	col = max_col;
	row = max_row;

	int weight = 0;
	int insertions = 0;
	int deletions = 0;

	while(row > 0 && col > 0)
	{
		int d = (rd1.seq[row] == rd2.seq[col]) ? MATCH : MISMATCH;
		int score = matrix[row][col];
		char arr = 'L';

		if(score == (matrix[row-1][col-1] + d))
			arr = 'D';
		else if(score == (matrix[row-1][col] + GAP))
			arr = 'U';

		switch(arr)
		{
			case 'D':

				if(rd1.seq[row] != rd2.seq[col] && rd1.seq[row] != 'N' && rd2.seq[col] != 'N')
					weight++;
				row--; col--;
				break;

			case 'L':

				weight++;
				deletions++;
				col--;
				break;

			case 'U':

				weight++;
				insertions++;
				row--;
				break;
		}
	}

	int ovl_size = int(((max_row - row) + (max_col - col))/2);
	int indels = insertions + deletions;

	isContained = false;	// assume it is false

	if( (!nogap || indels == 0) && ovl_size > min_overlap && ((double)weight/ovl_size) <= epsilon )
	{
		if(row == 0 && max_row == rd1.length-1) // a non-proper overlap
		{
			if(!(col == 0 && max_col == rd2.length-1 && rd1.id > rd2.id))	// if the reads are identical, spare the read with the larger id
				isContained = true;
		}

		theCache[numCached].IDs[0] = rd1.id;
		theCache[numCached].IDs[1] = rd2.id;
		theCache[numCached].strand = strand;
		theCache[numCached].weight = weight;
		theCache[numCached].indels = indels;

		if(strand)
		{
			theCache[numCached].pos[0] = rd1.length-max_row-1;
			theCache[numCached].pos[1] = rd2.length-max_col-1;
			theCache[numCached].pos[2] = rd1.length-row-1;
			theCache[numCached].pos[3] = rd2.length-col-1;
		}
		else
		{
			theCache[numCached].pos[0] = row;
			theCache[numCached].pos[1] = col;
			theCache[numCached].pos[2] = max_row;
			theCache[numCached].pos[3] = max_col;
		}

		theCache[numCached].score = ( (2.0 * theCache[numCached].weight) /
			(theCache[numCached].pos[2] - theCache[numCached].pos[0] + theCache[numCached].pos[3] - theCache[numCached].pos[1]) );

		numCached += 1;

		return true;
	}
	return false;
}
