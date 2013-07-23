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

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <climits>
#include <sstream>
#include <cassert>

#ifdef USEOPENMP
	#include <omp.h>
#endif

#include "KmerNode.h"
#include "KmerList.h"
#include "KmerIter.h"
#include "KmerHash.h"
#include "SmithWaterman.h"
#include "ReadMatch.h"
#include "Naive.h"
#include "DnaRead.h"
#include "ReadAln.h"
#include "HapUtils.h"

#include "Dna.h"

using namespace std;

void DNA::set_SmithWaterman(int max_rs, int min_overlap, int match_score, int mismatch_score,
	int indel_score, int gap_score, double epsilon, int max_quality, int num_nucleotides, int max_steps, double homozyg_prior)
{
	SW = new SmithWaterman(max_rs, min_overlap, match_score, mismatch_score,
		indel_score, epsilon, gap_score, max_quality, num_nucleotides, max_steps, homozyg_prior);
	NB = SW->NB;
}

int OCompare(const void * a, const void * b)
{
	if( (*(OvlCache*)a).score < (*(OvlCache*)b).score )
		return -1;
	else if( (*(OvlCache*)a).score > (*(OvlCache*)b).score )
		return 1;
	return 0;
}

int DNA::get_overlaps(int rd, int overlapLimit, ReadAln * rdalns, SmithWaterman * NW)
{
	int numOverlaps = NW->get_numCached();
	OvlCache * oc = NW->get_cache();

	if(numOverlaps > overlapLimit)
	{
		qsort(oc, NW->get_numCached(), sizeof(OvlCache), OCompare);
		numOverlaps = overlapLimit;
	}

	NW->NB->set_prior(0.999 - (0.001)*100*((numOverlaps+1.0)/(overlapLimit+1.0)));

	DnaRead read2(max_read_size);	// this is our chaperone read

	for(int j=0; j<numOverlaps; j++)
	{
		int rd2 = oc[j].IDs[1];
		int offset = oc[j].pos[0] - oc[j].pos[1];
		int width = oc[j].indels + 1;
		int next = 0;

		read2.replicate_with_quals(reads[rd2]);

		if(oc[j].strand)
			read2.revcomp(true);

		if(offset < 0)
		{
			rdalns[j].pos = NW->align(reads[rd], read2, 1, (0 - offset) - width, 2*width+1, next, rdalns[j], false);
			rdalns[j].start = (next >= 0) ? 0 : (0 - next);
		}
		else
		{
			rdalns[j].pos = NW->align(read2, reads[rd], 1, offset - width, 2*width+1, next, rdalns[j], true);
			rdalns[j].start = (next >= 0) ? next : 0;
		}
	}

	return numOverlaps;
}

bool DNA::handle_indels(int j, ReadAln * rdalns, int num_ovls, bool deletion)
{
	if(platform == ILLUMINA)
		return false;

	int yes = 0;
	int no = 1;		// the vote of the read itself

	for(int i=0; i<num_ovls; i++)
	{
		if(rdalns[i].start <= j && rdalns[i].pos >= 0 && rdalns[i].read1[ rdalns[i].pos ] != 'N') // the ambiguous letters do not get to vote!!
		{
			char ch = (deletion) ? rdalns[i].read1[ rdalns[i].pos ] : rdalns[i].read0[ rdalns[i].pos ];
			if(ch == '-')
				yes++;
			else
				no++;
		}
	}

	if( 3 > no && yes >= 3 )
		return true;

	return false;
}

void DNA::correct_read(ofstream & newReads, char * new_quals, int rd, ReadAln * rdalns, int num_ovls, double * votes)
{
	newReads << defs[rd] << "\n";

	double infi = 1000000.0;
	int bytes = 0;

	for(int j=0; j<reads[rd].length - 1; j++)	// j will keep track of the index of the original read
	{
		char orig = reads[rd].seq[j+1];					// the original base
		short int orig_score = ((short int) (reads[rd].quals[j+1])) - DnaRead::PHRED;	// the original score

		double min_error = infi;
		char winner = orig;

		int highest = 0;

		bool deletion = false;

		if( handle_indels(j, rdalns, num_ovls, false) )
		{
			min_error = infi;
			for(int k=0; k<4; k++)
			{
				char let = alph[k];
				votes[(int)let] = 0.0;		// reset the vote
				highest = 0;

				for(int i=0; i<num_ovls; i++)
				{
					if(rdalns[i].start <= j && rdalns[i].pos >= 0 && rdalns[i].read0[ rdalns[i].pos ] == '-')
					{
						if(rdalns[i].read1[ rdalns[i].pos ] == let)
						{
							votes[(int)let] += NB->log_same[ rdalns[i].qual1[rdalns[i].pos] ] [ rdalns[i].weight ];
							highest = (highest < rdalns[i].qual1[rdalns[i].pos]) ? rdalns[i].qual1[rdalns[i].pos] : highest;
						}
						else
							votes[(int)let] += NB->log_diff[ rdalns[i].qual1[rdalns[i].pos] ] [ rdalns[i].weight ];
					}
				}

				if(highest > 0 && votes[(int)let] < min_error)
				{
					min_error = votes[(int)let];
					winner = let;
				}
			}

			for(int i=0; i<num_ovls; i++)
			{
				if(rdalns[i].start <= j && rdalns[i].pos >= 0 && rdalns[i].read0[ rdalns[i].pos ] == '-')		// the other alignments are frozen
					rdalns[i].pos -= 1;
			}

			j--;			// since we will process this position again
		}
		else
		{
		    // if insertion didn't win, just remove all gaps up to the current base:
		    for(int i=0; i<num_ovls; i++)
            {
                if(rdalns[i].start <= j && rdalns[i].pos >= 0)
                {
                    while(rdalns[i].read0[ rdalns[i].pos ] == '-')
                        rdalns[i].pos -= 1;
                }
            }

			deletion = handle_indels(j, rdalns, num_ovls, true);
			if( !deletion )
			{
				min_error = infi;
				for(int k=0; k<4; k++)
				{
					highest = 0;
					char let = alph[k];
					if(orig == let)
					{
						highest = orig_score;
						votes[(int)let] = NB->log_correct[ orig_score ];		// this also resets the votes left over from previous round
					}
					else
						votes[(int)let] = NB->log_error[ orig_score ];

					for(int i=0; i<num_ovls; i++)
					{
						if(rdalns[i].start <= j && rdalns[i].pos >= 0)
						{
							assert( orig == rdalns[i].read0[ rdalns[i].pos ] );

							if(rdalns[i].read1[ rdalns[i].pos ] == let)
							{
								votes[(int)let] += NB->log_same[ rdalns[i].qual1[rdalns[i].pos] ] [ rdalns[i].weight ];
								highest = (highest < rdalns[i].qual1[rdalns[i].pos]) ? rdalns[i].qual1[rdalns[i].pos] : highest;
							}
							else if( rdalns[i].read1[ rdalns[i].pos ] == orig )
								votes[(int)let] += NB->log_error[ rdalns[i].qual1[rdalns[i].pos] ];
							else
								votes[(int)let] += NB->log_diff[ rdalns[i].qual1[rdalns[i].pos] ] [ rdalns[i].weight ];
						}
					}

					if(highest > 0 && ((votes[(int)let] < min_error) || (votes[(int)let] == min_error && let == orig)))
					{
						min_error = votes[(int)let];
						winner = let;
					}
				}
			}

			for(int i=0; i<num_ovls; i++)
			{
				if(rdalns[i].start <= j && rdalns[i].pos >= 0)
					rdalns[i].pos -= 1;
			}
		}

		if(!deletion)
		{
			newReads << winner;
			new_quals[bytes++] = (char) (orig_score + DnaRead::PHRED);
		}
	}

	newReads << "\n+\n";
	for(int i=0; i<bytes; i++)
		newReads << new_quals[i];
	newReads << "\n";
}

string DNA::read_next_fasta(ifstream & fh, DnaRead & rd, bool & more)
{
	string s;
	char ch;
	int pos = 1;
	more = false;

	while(getline(fh, s))
	{
		istringstream instr(s);
		instr >> ch;
		if(ch == '>')	// new read
		{
			more = true;
			break;
		}

		do
		{
			if(pos >= max_read_size)
				throw "Read is too long or fasta file is corrupted!";

			rd.seq[pos++] = capit[(int) ch ];
		}while(instr >> ch);
	}
	rd.set_length(pos);
	return s;
}

string DNA::read_next_fastq(ifstream & fh, DnaRead & rd, bool & more, bool get_quals)
{
	string s;
	char ch;
	int pos = 1;

	getline(fh, s);
	if((int)s.size() >= max_read_size)
		throw "Read is too long or fastq file is corrupted!";

	istringstream instr(s);

	while(instr >> ch)
		rd.seq[pos++] = capit[(int) ch ];

	rd.set_length(pos);

	getline(fh, s);	// for the second defline
	getline(fh, s); // for the qualities

	if((int)s.size() >= max_read_size)
		throw "Number of quality values does not match the read length in fastq file!";

	if(get_quals)
	{
		istringstream instr2(s);
		pos = 1;
		while(instr2 >> ch)
		{
			if(((int)ch - DnaRead::PHRED) < 0) throw "Negative quality values are not allowed! Check PHRED_OFFSET.";
			rd.quals[pos++] = ch;
		}

		if(pos != rd.length)
			throw "Number of quality values does not match the read length in fastq file!";
	}

	more = false;
	if(getline(fh, s))
		more = true;	// there are more reads

	return s;
}

void DNA::prune_overlaps(char * filename, int nthreads)
{
	int rd1, rd2, strand, st1, st2, end1, end2, diff, indels;

	bool * dups = new bool[num_reads+1];
	bool * starts = new bool[num_reads+1];
	bool * ends = new bool[num_reads+1];
	for(int i=0; i<=num_reads; i++)
	{
		dups[i] = false;
		starts[i] = false;
		ends[i] = false;
	}

	for(int tid=0; tid<nthreads; tid++)
	{
		int remainder = tid%alphabet_size;
		int dividend = (tid-remainder)/alphabet_size;

		string s(filename);
		s += ".ovltmp";
		s += (char)('A' + dividend);
		s += (char)('A' + remainder);

		char * fname = new char[s.size()+1];
		strcpy(fname, s.c_str());

		ifstream ovlFile;
		open_n_check(ovlFile, fname);

		while(ovlFile >> rd1 >> rd2 >> strand >> st1 >> st2 >> end1 >> end2 >> diff >> indels)
		{
			if(st1==0 && st2==0 && end1==reads[rd1].length-1 && end2==reads[rd2].length-1)	// duplicates
			{
				if(rd1 < rd2)
					dups[rd1] = true;	// for duplicates always remove the read with the smaller id
				else
					dups[rd2] = true;
			}
			else if(st1==0 && end1==reads[rd1].length-1)	// rd1 is contained
				dups[rd1] = true;
			else if(st2==0 && end2==reads[rd2].length-1)
				dups[rd2] = true;
			else
			{									// note that if we are here the overlap must be proper!
				if(st1 == 0)
					starts[rd1] = true;
				else
					ends[rd1] = true;

				if(strand == 0)
				{
					if(st2 == 0)
						starts[rd2] = true;
					else
						ends[rd2] = true;
				}
				else
				{
					if(st2 == 0)
						ends[rd2] = true;
					else
						starts[rd2] = true;
				}
			}
		}
		check_n_close(ovlFile);
		delete [] fname;
	}

	for(int i=0; i<=num_reads; i++)
	{
		if(!starts[i] || !ends[i])	// if the read is an end node, mark it for deletion
			dups[i] = true;
	}

	string s2(filename);
	s2 += ".overlaps";

	char * fname2 = new char[s2.size()+1];
	strcpy(fname2, s2.c_str());

	ofstream after;
	open_n_check(after, fname2);

	for(int tid=0; tid<nthreads; tid++)
	{
		int remainder = tid%alphabet_size;
		int dividend = (tid-remainder)/alphabet_size;

		string s(filename);
		s += ".ovltmp";
		s += (char)('A' + dividend);
		s += (char)('A' + remainder);

		char * fname = new char[s.size()+1];
		strcpy(fname, s.c_str());

		ifstream before;
		open_n_check(before, fname);

		while(before >> rd1 >> rd2 >> strand >> st1 >> st2 >> end1 >> end2 >> diff >> indels)
		{
			if( rd1 < rd2 && !dups[rd1] && !dups[rd2] )
				after<<rd1<<" "<<rd2<<" "<<strand<<" "<<st1<<" "<<st2<<" "<<end1<<" "<<end2<<" "<<diff<<" "<<indels<< "\n";
		}
		check_n_close(before);
		delete [] fname;
	}
	check_n_close(after);

	delete [] fname2;
	delete [] dups;
	delete [] ends;
	delete [] starts;
}

void DNA::cat_reads(char * filename, int nthreads)
{
	if(nthreads == 1)
		return;

	ofstream outFile;
	open_n_check(outFile, filename);

	int buffsize = 32*1024*1024;
	char * buff = new char[buffsize];

	for(int tid=0; tid<nthreads; tid++)
	{
		int remainder = tid%alphabet_size;
		int dividend = (tid-remainder)/alphabet_size;

		string s2(filename);
		s2 += ".cortmp";
		s2 += (char)('A' + dividend);
		s2 += (char)('A' + remainder);

		char * fname = new char[s2.size()+1];
		strcpy(fname, s2.c_str());

		ifstream inFile;
		open_n_check(inFile, fname);

		while(true)
		{
			inFile.read(buff, buffsize);
			int c = inFile.gcount();
			outFile.write(buff, c);
			if(!inFile) break;
		}
		check_n_close(inFile);
	}

	delete [] buff;
	check_n_close(outFile);
}

int DNA::write_contig(ofstream & hapFile, ofstream & mapFile, ContigTuple * contups, int contig_num, int tuple_num)
{
	hapFile << ">contig_" << contig_num << "\n";

	int rd1 = contups[0].read;

	if(contups[0].strand == 1)
		reads[rd1].revcomp(true);
	else
		reads[rd1].revcomp(false);

	ReadAln * rdalns = new ReadAln[tuple_num];
	for(int i=0; i<tuple_num; i++)
		rdalns[i].set(2*SW->max_read_size);

	int count = -1;
	int start = -1;
	int textwidth = 80;

	int bytes = 0;
	double * votes = new double[256];
	int i=0, j=1, rd2=0;

	while( i < tuple_num - 1 )
	{
		j = i+1;

		rd1 = contups[i++].read;
		rd2 = contups[j].read;

		mapFile << rd1 <<" "<< contups[i-1].strand <<" "<< contig_num <<" "<< bytes+1 <<" "<< reads[rd1].length-1 << " -1\n";

		if(contups[j].strand == 1)
			reads[rd2].revcomp(true);
		else
			reads[rd2].revcomp(false);

		int next = 0;
		count++;
		int width = (reads[rd1].length > reads[rd2].length) ? (int) (reads[rd1].length * SW->epsilon) : (int) (reads[rd2].length * SW->epsilon);
		rdalns[count].pos = SW->align(reads[rd2], reads[rd1], 1, contups[j].offset - width, 2*width + 1, next, rdalns[count], false, true);

		if(rdalns[count].pos <= 0 || next < 0)
			rdalns[count].pos = SW->align(reads[rd2], reads[rd1], 1, 0-reads[rd1].length, 2 * reads[rd1].length , next, rdalns[count], false, true);

		assert(next >= 0);
		assert(rdalns[count].pos > 0);

		if(count == 0)
			reads[rd1].write_seq(hapFile, 1, next, bytes, textwidth);

		if(count - start > 1)
		{
			char * letters = new char[2 * (count - start)];
			short int * scores = new short int[2 * (count - start)];

			for(int n=0; n<2*(count-start); n++)				// initialize with gaps
			{
				letters[n] = '-';
				scores[n] = 0;
			}

			int letpos = 0;

			for(int ind=0; ind<next; ind++)		// write the consensus till the next offset
			{
				letpos = 0;

				int p = rdalns[count-1].pos--;
				assert(p >= 0);
				char let1 = rdalns[count-1].read1[p];
				char let0 = rdalns[count-1].read0[p];

				short int score1;
				short int score0;

				scores[letpos] = rdalns[count-1].qual1[p];
				letters[letpos++] = let0;

				scores[letpos] = rdalns[count-1].qual0[p];
				letters[letpos++] = let1;

				if( let0 == '-' )		// a gap in the last read, just increase the offset
					assert(let1 != '-');

				for(int k=count-2; k>start; k--)	// note that we traverse from the last couple to the first!!
				{
					p = rdalns[k].pos;
					assert(p >= 0);

					if(let1 == '-')		// if the current read had a gap in the previous alignment
					{
						let1 = rdalns[k].read1[p];
						let0 = rdalns[k].read0[p];

						score1 = rdalns[k].qual1[p];
						score0 = rdalns[k].qual0[p];

						if(let0 == '-')		// if it's current position is also a gap there is nothing to do, just proceed
						{
							rdalns[k].pos -= 1;
						}
						else	// otherwise we need to "insert a gap" into this alignment (we pretend to do so by keeping the cursor where it was)
						{
							let0 = '-';
							let1 = '-';
							score0 = SW->gap_score;
							score1 = SW->gap_score;
						}
					}
					else	// if it wasn't a gap in the previous alignment
					{
						let1 = rdalns[k].read1[p];
						let0 = rdalns[k].read0[p];

						score1 = rdalns[k].qual1[p];
						score0 = rdalns[k].qual0[p];

						if(let0 == '-')	// but it's a gap now, that means everything we processed so far must also be gaps
						{
							int temp = letpos-1;
							for(int m=k+1; m<count; m++)
							{
								if(letters[temp] == '-' && letters[temp-1] == '-')
								{
									break;
								}
								else
								{
									rdalns[m].pos++;	// else step back
								}

								scores[temp] = SW->gap_score;
								letters[temp--] = '-';

								scores[temp] = SW->gap_score;
								letters[temp--] = '-';
							}
						}
						rdalns[k].pos -= 1;
					}

					scores[letpos] = score0;
					letters[letpos++] = let0;

					scores[letpos] = score1;
					letters[letpos++] = let1;
				}

				votes[(int)'A'] = 0;
				votes[(int)'C'] = 0;
				votes[(int)'G'] = 0;
				votes[(int)'T'] = 0;
				votes[(int)'-'] = 0;

				// then we can get the consensus
				int num_of_votes = 0;
				votes[(int) letters[0] ] += scores[0];
				num_of_votes += scores[0];

				if(letters[0] == '-')
					ind--;

				for(int n=1; n<letpos; n+=2)
				{
					votes[(int) letters[n] ] += scores[n];

					if(n+1 < letpos)
						assert(letters[n] == letters[n+1]);

					num_of_votes += scores[n];
				}

				char best = '-';
				double max_vote = 0.0;
				for(int z=0; z<=4; z++)
				{
					char let = alph[z];
					if(votes[(int)let] > max_vote)
					{
						max_vote = votes[(int)let];
						best = let;
					}
				}

				int inc = 1;

				if(best != '-')
					hapFile << best;
				else
					inc = 0;

				bytes += inc;
				if(inc && bytes%textwidth == 0)
					hapFile << "\n";

				for(int k=start+1; k<count; k++)
				{
					if(rdalns[k].pos < 0)
					{
						start = k;
						if(count - start < 2)	// it means the current region is only covered by a single read
						{
							// we are at position (ind) and we have to output (next - ind) more letters
							for(int s=ind+2; s<=next; s++)
							{
								assert(reads[rd1].length >= s);
								hapFile << reads[rd1].seq[s];	// note that rd2 is not among the alignments yet!!

								bytes++;
								if(bytes%textwidth == 0)
									hapFile << "\n";
							}
							next = -1;	// to break out the outer for loop
							break;
						}
					}
				}
			}
			delete [] letters;
			delete [] scores;
		}
	}

	if(rd2>0)
	{
		mapFile << rd2 <<" "<< contups[j].strand <<" "<< contig_num <<" "<< bytes+1 <<" "<< reads[rd2].length-1 << " -1\n";
		reads[rd2].write_seq(hapFile, 1, reads[rd2].length-1, bytes, textwidth);
	}

	delete [] rdalns;
	delete [] votes;
	return bytes;
}

bool * DNA::scan_contigs(char * inputFilename, int minUnique)
{
	int offset, strand, read, dummy;
	int * dups = new int[num_reads+1];
	bool * conts = new bool[num_reads+1];
	for(int i=0; i<=num_reads; i++)
	{
		dups[i] = 0;
		conts[i] = false;
	}

    ifstream contigFile;
	open_n_check(contigFile, inputFilename);

	while(contigFile >> offset)		// read the contig file
	{
		if(offset == -1)	// then it means we have a new contig
			contigFile >> dummy;		// reserved
		else if(offset == -2)
		{
			int rd, same, begin, end;
			contigFile >> rd >> same >> begin >> end;	// ignore contained reads for now
		}
		else
		{
			contigFile >> strand >> read;
			dups[read] += 1;
		}
	}

	int contigNum = 1, duplicates = 0, total = 0;
	check_n_close(contigFile);

	ifstream contigFile2;
	open_n_check(contigFile2, inputFilename);
	read = 0;
	while(contigFile2 >> offset)		// read the contig file
	{
		if(offset == -1)	// then it means we have a new contig
		{
			contigFile2 >> dummy;		// reserved
			conts[contigNum++] = (total - duplicates < minUnique) ? false : true;
			total = 0;
			duplicates = 0;
			read = 0;
		}
		else if(offset == -2)
		{
			int rd, same, begin, end;
			contigFile2 >> rd >> same >> begin >> end;	// ignore contained reads for now
		}
		else
		{
			if(dups[read] > 1)		// this is the previous read unless it's the first
				duplicates += offset;
			contigFile2 >> strand >> read;
			total += offset;
		}
	}

	check_n_close(contigFile2);
	return conts;
}

void DNA::write_contigs(char * inputFilename, int minContigSize, char * outputFilename)
{
    bool * conts = scan_contigs(inputFilename, minContigSize);

    string s(outputFilename);		// this is an auxiliary file, that can be used for analysis or as input to a scaffolder
	s += ".tmp";
	char * fname = new char[s.size()+1];
	strcpy(fname, s.c_str());
	ofstream mapFile;
	open_n_check(mapFile, fname);

	ofstream hapFile;
	open_n_check(hapFile, outputFilename);

	ifstream contigFile;
	open_n_check(contigFile, inputFilename);

	int offset, strand, read;
	int contigNum = 1;
	ContigTuple * contups = new ContigTuple[num_reads];
	int * contig_len = new int[num_reads];
	int tuple_num = 0;
	int len = 0;
	long int totalLength = 0;
	int dummy;

	int numOfContigs = 0;

	while(contigFile >> offset)		// read the contig file
	{
		if(offset == -1)	// then it means we have a new contig
		{
			contigFile >> dummy;		// reserved
			if(dummy == 0 && tuple_num > 1 && conts[contigNum] && (len + reads[contups[tuple_num-1].read].length - 1) > minContigSize)
			{
				// returns the total number of bytes written to the file (for one haplome only)
				len = write_contig(hapFile, mapFile, contups, contigNum, tuple_num);

				if(len%80!=0)
					hapFile << "\n";

				totalLength += len;
				contig_len[numOfContigs++] = len;

				if(numOfContigs == num_reads)
					throw "Contig number exceeds number of reads!";
			}

			tuple_num = 0;
			len = 0;
			contigNum++;
		}
		else if(offset == -2)
		{
			int rd, same, begin, end;
			contigFile >> rd >> same >> begin >> end;	// ignore contained reads for now
		}
		else
		{
			contigFile >> strand >> read;
			if(tuple_num > 0)	// do not count the first offset
				len += offset;
			if(tuple_num >= num_reads)
				throw "Reads are duplicated in contig! Aborting...";

			contups[tuple_num].offset = offset;
			contups[tuple_num].read = read;
			contups[tuple_num++].strand = (strand & (OUTIE|INNIE));
		}
	}
	check_n_close(contigFile);

	length_stats(contig_len, numOfContigs, totalLength);

	delete [] contups;
	delete [] contig_len;
	delete [] conts;
	delete [] fname;

	check_n_close(hapFile);
	check_n_close(mapFile);
}

DNA::DNA(int kmer_size, int max_read_size, int phred_offset, Sequencer platform)
{
	SW = NULL;
	NB = NULL;
	KH = NULL;

	this->max_read_size = max_read_size;
	this->kmer_size = kmer_size;
	this->platform = platform;

	reads = NULL;
	deflines = NULL;
	sequences = NULL;
	defs = NULL;

	num_reads = 0;

	// A = 00, C = 01, G = 10, T = 11
	alphabet = new int[256];	// this will serve as a hash
	alphabet[(int)'A'] = 0;
	alphabet[(int)'C'] = 1;
	alphabet[(int)'G'] = 2;		// note that we require all the bases to be capital, this has to be ensured when the reads are being read!!
	alphabet[(int)'T'] = 3;
	alphabet[(int)'N'] = 4;

	capit = new char[256];	// this will serve to capitalize reads
	for(int i=0; i<256; i++)
		capit[(int)i] = 'N';		// this ensures any letter other than A,C,G,T is converted to N

	capit[(int)'A'] = 'A';
	capit[(int)'C'] = 'C';
	capit[(int)'G'] = 'G';
	capit[(int)'T'] = 'T';
	capit[(int)'a'] = 'A';
	capit[(int)'c'] = 'C';
	capit[(int)'g'] = 'G';
	capit[(int)'t'] = 'T';

	alph = new char[6];
	alph[0] = 'A';
	alph[1] = 'C';
	alph[2] = 'G';
	alph[3] = 'T';
	alph[4] = '-';
	alph[5] = '+';

	DnaRead::set_phred(phred_offset);
	expected_overlaps = 1;
	rmax_overlaps = 1;

	isConsensr = false;
	isCorrectr = false;
	isOverlappr = false;
}

DNA::~DNA()
{
	if(reads != NULL)
	{
		for(int i=0; i<=num_reads; i++)
			reads[i].detach();

		delete [] reads;
		delete [] sequences;
	}

	if(defs != NULL)
	{
		delete [] defs;
		delete [] deflines;
	}
	if(KH != NULL)
		delete KH;

	reads = NULL;
	defs = NULL;
	KH = NULL;

	delete [] alphabet;
	delete [] capit;
	delete [] alph;
	delete SW;
}

void DNA::trim_reads(char * reads_filename, char * reads_filename2, char * output_filename,
	int winsize, double threshold, int revcomp)
{
	string s1;
	ifstream readsFile;
	open_n_check(readsFile, reads_filename);
	getline(readsFile, s1);

	string s2;
	ifstream readsFile2;
	bool paired = false;
	if(reads_filename2 != NULL)
	{
		open_n_check(readsFile2, reads_filename2);
		getline(readsFile2, s2);
		paired = true;
	}

	ofstream outFile;
	bool useStdout = true;
	if(output_filename != NULL)
	{
		open_n_check(outFile, output_filename);
		useStdout = false;
	}

	DnaRead rd(max_read_size);

	int max_trimmed_size = 0;
	long int after = 0;
	long int before = 0;

	bool trim = true;

	if(threshold > 0.0001)
		threshold = 0 - log10(threshold);
	else
		trim = false;

	int read_id = 1;
	bool more = true;
	bool more2 = false;
	int order = 1;

	string s;
	string s1new;
	string s2new;

	int droppoint = 1+2*winsize;

	while( more || more2 )
	{
		if(!paired || order == 1)
		{
			s1new = read_next_fastq(readsFile, rd, more, true);
			s = s1;
			s1 = s1new;
		}
		else
		{
			s2new = read_next_fastq(readsFile2, rd, more2, true);
			s = s2;
			s2 =s2new;
		}

		rd.id = read_id++;

		if(rd.length >= max_read_size)
			throw "Read length exceeds maximum!";

		before += (rd.length-1);

		if( trim && rd.trim(winsize, droppoint, NB, threshold) )
			after += (rd.length-1);

		bool opposite = (revcomp & order) ? true : false;

		if(useStdout)
			rd.write_fastq(s, opposite);
		else
			rd.write_fastq(outFile, s, opposite);

		if(rd.length > max_trimmed_size)
			max_trimmed_size = rd.length;

		order = (order == 1) ? 2 : 1;
		rd.rev = false;	  // reset the state
	}

	if(!trim) after = before;

	cerr << "Number of reads : " << read_id - 1 << endl;
	cerr << "Maximum trimmed read length : " << max_trimmed_size-1 << endl;
	cerr << "Total number of bases before trimming : " << before << endl;
	cerr << "Total number of bases after trimming : " << after << endl;

	check_n_close(readsFile);
	if(paired) check_n_close(readsFile2);
	if(!useStdout) check_n_close(outFile);
}

int DNA::read_reads(char * reads_filename, long int genome, bool getQuals, bool getDeflines, char * output_filename)
{
	long int num = 0;
	long int sum = 0;
	long int defsum = 0;

	int buffsize = 32*1024*1024;
	char * buff = new char[buffsize];

	ifstream inFile;
	open_n_check(inFile, reads_filename);
	int ctr = 0;
	long int linelen = 0;
	int max_rs = 0;

	while(true)
	{
		inFile.read(buff, buffsize);
		int c = inFile.gcount();
		for(int i=0; i<c; i++)
		{
			if(buff[i] == '\n')
			{
				ctr++;
				if(ctr==1)
				{
					defsum += linelen;
				}
				if(ctr==2)	// read the sequence
				{
					num++;
					sum += linelen;
					if(linelen > max_rs)
						max_rs = linelen;
				}
				else if(ctr==4)
					ctr = 0;

				linelen = 0;
			}
			else
				linelen += 1;
		}

		if(!inFile) break;
	}
	check_n_close(inFile);

	max_read_size = max_rs + 5;
	SW->set_matrix(max_read_size);

	if(num > (2147483600/2)) throw "Number of reads exceeds the limit of one billion. Program is aborting";
	if(ctr != 0) throw "Number of lines in fastq file is not a multiple of 4. File may be corrupted!";
	if(sum == 0) throw "Fastq file is empty or corrupted!";

	string s, s2;

	num_reads = (int)num;
	double coverage = 0.0;

	if(!isConsensr)
	{
		if(genome > 270000000)
			this->kmer_size = 15;

		coverage = (sum-num_reads)/(double)(genome);
		double average = (sum-num_reads)/(double)num_reads;	// average read length
		double minOvl = average - (average * log(genome))*((double)genome /(sum-num_reads));
		minOvl = 0.95*minOvl;
		if(minOvl < 0)
		{
			SW->min_overlap = 1+2*this->kmer_size;
			cerr << "WARNING: Coverage is insufficient!" << endl;
		}
		else if( minOvl < 1+2*kmer_size)
		{
			this->kmer_size = max(11, int(minOvl/ 2));		// do not allow a kmer-size less than 11
			SW->min_overlap = 1+2*this->kmer_size;
		}
		else	// all is fine
		{
			SW->min_overlap = int(minOvl);
		}
		cerr << "Setting the minimum overlap to : " << SW->min_overlap << endl;
		cerr << "Setting the kmer size to : " << this->kmer_size << endl;
	}

	sequences = (getQuals) ? new char[2*(sum+num_reads+1)] : new char[sum+num_reads+1];
	if(getDeflines)
	{
		deflines = new char[defsum+num_reads+2];
		defs = new char*[num_reads+2];
		defs[1] = &deflines[0];
	}

	reads = new DnaRead[num_reads+1];
	DnaRead tmp(max_read_size);

	ofstream outputFile;
	if(output_filename != NULL)
	{
		string s(output_filename);
		s += ".reads";
		char * fname = new char[s.size()+1];
		strcpy(fname, s.c_str());

		open_n_check(outputFile, fname);
		outputFile << num_reads << "\n";
		delete [] fname;
	}

	ifstream readsFile;
	open_n_check(readsFile, reads_filename);

	linelen = 1;
	ctr = 0;

	long int pos = 0;
	long int defpos = 0;
	int n = 1;
	reads[n].seq = &sequences[pos];
	reads[n].seq[0] = '-';
	reads[n].id = n;

	while(true)
	{
		readsFile.read(buff, buffsize);
		int c = readsFile.gcount();
		for(int i=0; i<c; i++)
		{
			if(buff[i] == '\n')
			{
				ctr++;
				if(ctr==1 && getDeflines)	// read the defline
				{
					deflines[defpos++] = '\0';
					defs[n+1] = &deflines[defpos];
				}

				if(ctr==2)					// read the sequence
				{
					pos += linelen;
					reads[n].length = linelen;

					if(output_filename != NULL)
						outputFile << linelen - 1 << "\n";

					if(getQuals)
					{
						reads[n].quals = &sequences[pos];
						reads[n].quals[0] = (char) (DnaRead::PHRED + 1);
					}
				}
				else if(ctr==4)				// read the qualities (and done with the current read)
				{
					if(getQuals)
					{
						pos += linelen;
						assert(reads[n].length == linelen);
					}

					if(n < num_reads)
					{
						reads[++n].seq = &sequences[pos];		// assing the memory for the next read
						reads[n].seq[0] = '-';					// init the empty base
						reads[n].id = n;
					}
					ctr = 0;				// end of the cycle
				}
				linelen = 1;
			}
			else if(ctr==0 && getDeflines)
				deflines[defpos++] = buff[i];
			else if(ctr==1)
				reads[n].seq[linelen++] = capit[(int)buff[i]];
			else if(ctr==3 && getQuals)
				reads[n].quals[linelen++] = buff[i];
		}
		if(!readsFile) break;
	}

	if(output_filename != NULL)
		check_n_close(outputFile);

	check_n_close(readsFile);
	delete [] buff;

	return int(coverage);
}


void DNA::fill_in_kmers(long int total_kmers, int max_kmer_duplicity, DnaRead * dnareads, int num)
{
	long int hash_size = 1 << (2*kmer_size);
	long int max_nodes = total_kmers;

	KH = new KmerHash(hash_size, max_nodes, max_kmer_duplicity);

	for(int rd=1; rd<=num; rd++)
	{
		int len = dnareads[rd].length;

		for(int i=1; i<len; i+=kmer_size)
		{
			int kmer_value = 0;
			bool ok = true;

			if(i+kmer_size > len)
				i = len-kmer_size;

			for(int j=0; j<kmer_size; j++)
			{
				int ch = alphabet[(int) dnareads[rd].seq[i+j] ];
				if(ch > 3)
				{
					ok = false;		// do not allow kmers with N
					break;
				}
				kmer_value = ((kmer_value << 2)  + ch);
			}
			if(ok)
			{
				KH->add_entry(kmer_value, rd, (short int) i);
			}
		}
	}
}

int RCompare(const void * a, const void * b)
{
	int c = (*(RMatch*)a).ID - (*(RMatch*)b).ID;
	if(c != 0)
		return c;

	return (*(RMatch*)a).hit - (*(RMatch*)b).hit;
}


bool DNA::compare_reads(DnaRead & read1, bool strand, ReadMatch * RM, RMatch * rd_list, DnaRead * dnareads, SmithWaterman * NW)
{
	int len = read1.length;
	int num_pairs = 0;
	int kmer_values[4];

	for(int i=1; i+kmer_size<=len; i++)
	{
		kmer_values[0] = 0;
		int numk = 1;

		bool ok = true;
		for(int j=0; j<kmer_size; j++)
		{
			int ch = alphabet[(int) read1.seq[i+j] ];

			if(ch > 3)		// the kmer contains an N
			{
				if(numk==1)		// then branch out
				{
					numk = 4;
					for(int k=numk-1; k>=0; k--)	// so that kmer_values[0] is changed last
					{
						ch = alphabet[(int) alph[k] ];
						kmer_values[k] = (kmer_values[0] << 2) + ch;
					}
				}
				else	// already branched out, kill the kmer
				{
					ok = false;
					break;
				}
			}
			else
			{
				for(int k=0; k<numk; k++)
					kmer_values[k] = (kmer_values[k] << 2) + ch;
			}
		}

		if(ok)
		{
			for(int k=0; k<numk; k++)
			{
				if( KH->kmer_list[kmer_values[k]].get_size() < KH->max_size )
				{
					KmerIter iter( &(KH->kmer_list[kmer_values[k]]) );
					for(KmerNode * kn = iter.get_first(); kn != NULL; kn = iter.get_next())
					{
						int rd2 = kn->data;
						short int diff = i - kn->pos;

						if( dnareads[rd2].id != read1.id )
						{
							rd_list[num_pairs].ID = dnareads[rd2].id;
							rd_list[num_pairs].count = 1;
							rd_list[num_pairs++].hit = diff;
						}
					}
				}
			}
		}
	}

	if(num_pairs == 0)
		return false;

	double eps = NW->epsilon;
	bool nogap = false;

	if(platform == ILLUMINA)
	{
		eps = 0.0;
		nogap = true;
	}

	qsort(rd_list, num_pairs, sizeof(RMatch), RCompare);

	int prevRead = rd_list[0].ID;
	int prevHit = rd_list[0].hit;
	int hitCountr = 1;
	int posCountr = 0;

	for(int i=1; i<num_pairs; i++)
	{
		if(rd_list[i].ID == prevRead && rd_list[i].hit == prevHit)
			hitCountr++;
		else
		{
			rd_list[posCountr].ID = prevRead;
			rd_list[posCountr].hit = prevHit;
			rd_list[posCountr++].count = hitCountr;

			prevHit = rd_list[i].hit;
			prevRead = rd_list[i].ID;
			hitCountr = 1;
		}
	}

	rd_list[posCountr].ID = prevRead;
	rd_list[posCountr].hit = prevHit;
	rd_list[posCountr++].count = hitCountr;

	int rd2 = rd_list[0].ID;
	for(int i=0; i<posCountr; i++)
	{
		if(rd_list[i].ID != rd2)
		{
			if(align_reads(read1, dnareads[rd2], strand, nogap, eps, RM, NW))
				return true;

			rd2 = rd_list[i].ID;
		}

		if(RM->num_hits < RM->max_hits)
		{
			RM->hits[RM->num_hits] = rd_list[i].hit;
			RM->counts[RM->num_hits++] = rd_list[i].count;
		}
	}

	return align_reads(read1, dnareads[rd2], strand, nogap, eps, RM, NW);
}

bool DNA::align_reads(DnaRead & read1, DnaRead & read2, bool strand, bool nogap, double eps, ReadMatch * RM, SmithWaterman * NW)
{
	bool isContained = false;

	if(RM->num_hits >= RM->max_hits)
		NW->align(read1, read2, 1, 0-read2.length, 2 * read2.length, strand, isContained, nogap);
	else
	{
		int k = 0;
		int k1, k2, leeway;

		RM->compress(NW->epsilon, kmer_size, read1.length, read2.length, nogap);	// this will bundle close hits into a single hit point

		while( RM->get_next_pos(k, eps, k1, k2, read1.length, read2.length, leeway, NW->min_overlap) )
		{
			if( NW->align(read1, read2, k1, k2, leeway, strand, isContained, nogap) ) // note: leeway is always at least 3
				break;
		}
	}
	RM->reset();	// clean up the hit list

	if(isContained && isOverlappr)	// return true only if running overlappr
		return true;

	return false;
}

int DNA::overlap_reads(int nthreads, char * filename, bool onestrand, int coverage)
{
	coverage = max(10, coverage);	// do not allow very small coverage values
	int coverage_cutoff = (coverage*100)/kmer_size;

	long int total_kmers = 0;
	for(int i=1; i<=num_reads; i++)
		total_kmers += (int(reads[i].length / (double)kmer_size) + 2);

	fill_in_kmers(total_kmers, coverage_cutoff, reads, num_reads);

	int max_overlaps = coverage_cutoff * max_read_size;		// it is not possible to find more pairs per read
	if(max_overlaps > num_reads) max_overlaps = num_reads;

	int tid;

#ifdef USEOPENMP

	omp_set_num_threads(nthreads);
	int step = (num_reads / nthreads) + 1;

	#pragma omp parallel private(tid)
	{
		tid = omp_get_thread_num();

		int st = 1 + (tid * step);
		int end = (st + step - 1 < num_reads) ? st + step - 1 : num_reads;
#else
		int st = 1;
		int end = num_reads;
		tid = 0;
		nthreads = 1;
#endif

		int remainder = tid%alphabet_size;
		int dividend = (tid-remainder)/alphabet_size;

		ofstream outputFile;

		char * newQuals = NULL;
		ReadAln * rdalns = NULL;
		double * votes = NULL;
		int overlapLimit = 1;

		if(isOverlappr)
		{
			string s(filename);
			s += ".ovltmp";
			s += (char)('A' + dividend);
			s += (char)('A' + remainder);

			char * fname = new char[s.size()+1];
			strcpy(fname, s.c_str());
			open_n_check(outputFile, fname);

			delete [] fname;
		}
		else
		{
			if(nthreads == 1)
				open_n_check(outputFile, filename);
			else
			{
				string s2(filename);
				s2 += ".cortmp";
				s2 += (char)('A' + dividend);
				s2 += (char)('A' + remainder);

				char * ofname = new char[s2.size()+1];
				strcpy(ofname, s2.c_str());
				open_n_check(outputFile, ofname);

				delete [] ofname;
			}

			newQuals = new char[2*max_read_size];
			overlapLimit = ((100 * coverage) > num_reads) ? num_reads : 100 * coverage;

			rdalns = new ReadAln[overlapLimit];
			for(int i=0; i<overlapLimit; i++)
				rdalns[i].set(2*max_read_size);

			votes = new double[256];
		}

		RMatch * rd_list = new RMatch[num_reads];
		ReadMatch * RM = new ReadMatch(14);

		DnaRead rd(max_read_size);
		SmithWaterman * NW = new SmithWaterman(*SW);
		NW->set_cache(max_overlaps);

		for(int i=st; i<=end; i++)
		{
			rd.replicate(reads[i]);
			bool isContained = compare_reads(rd, false, RM, rd_list, reads, NW);
			if(!onestrand && !isContained)
			{
				rd.revcomp();
				isContained = compare_reads(rd, true, RM, rd_list, reads, NW);
			}

			if(isOverlappr)
			{
				if(isContained)
					NW->reset_cache();
				NW->flush_cache(outputFile);
			}
			else
			{
				int counter = get_overlaps(i, overlapLimit, rdalns, NW);
				correct_read(outputFile, newQuals, i, rdalns, counter, votes);
				NW->reset_cache();
			}
		}

		delete NW;
		delete RM;
		delete [] rd_list;

		check_n_close(outputFile);

		if(isCorrectr)
		{
			delete [] rdalns;
			delete [] newQuals;
			delete [] votes;
		}

#ifdef USEOPENMP
	}
#endif

	return nthreads;
}
