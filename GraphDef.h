/*****************************************************************************
    $Author: Nilgun Donmez $
    $Date: 2012-09-20 17:27:54 -0400 (Thu, 20 Sep 2012) $

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

#ifndef GRAPH_DEF_H
#define GRAPH_DEF_H

#include <stdint.h>
#include <cstdlib>

// generic flag
#define NONE 0

#ifndef NULL
	#define NULL 0
#endif

// flags for arrows (up to 128)
#define OUTIE 1		// outward arrow
#define INNIE 2		// inward arrow
#define NON_PROPER 4	// a non-proper edge
#define FAKE 8
#define SUPPORTED 32
#define OK 64

// flags for nodes (up to 32000)
#define CONTAINED 1
#define DEAD 2
#define VISITED 4
#define REMOVED 256
#define LOOPY 128
#define REMOVABLE 16		// also used for edge flags
#define FOUND 512
#define SINGLE 1024

enum Sequencer { ILLUMINA, FOURFIVEFOUR, NOPLATFORM };
enum ClassID { READ_NODE = -1, MATE_PAIR = -2, PATH_NODE = -3, PATH_EDGE = -4 };

const int INDEG = 1;
const int OUTDEG = 2;
const int DEG = 3;

//**************************************************************
//	Typedefs
//**************************************************************

typedef uint8_t EdgeFlag;
typedef uint8_t ArrowType;
typedef short Direction;

typedef short NodeFlag;
typedef int NodeLength;

typedef int PathLength;
typedef short EdgeLength;

typedef char hapchar;

//**************************************************************
//	Structs
//**************************************************************

struct SDNode
{
	int target;
	EdgeLength in;
	EdgeLength out;
};

struct QStruct
{
	int readID;
	Direction dir;
	EdgeLength val;
};

struct PStruct
{
	int x[3];
};

struct ContigTuple
{
	short int offset;
	short int strand;
	int read;
};

struct OvlTag
{
	int rd2;
	int offset;
	int width;
	int strand;
	float err;
};

struct ReadTag
{
	int contigID;
	short int strand;
	int offset;
	int weight;
};

struct ContigTag
{
	int length;
	long int count;
};

struct OvlCache
{
	bool strand;
	float score;
	int IDs[2];
	EdgeLength pos[4];
	EdgeLength weight;
	EdgeLength indels;
};

struct RMatch
{
	int ID;
	short int hit;
	short int count;
};

struct ReadPair
{
	short int length;
	short int lib;
};

struct Library
{
	int insertSize;
	short int deviation;
	double sdwidth;
	int maxDist;
	int orient;

	double coverage;
	int minSupport;
};

//***************************************************************
//	Inline functions
//***************************************************************

inline int rev(int num) { return ( (num+1)%2 ); };

#endif // GRAPH_DEF_H
