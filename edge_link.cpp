#include "edge_link.h"
#include <stdlib.h>
#include <iostream>

edge_link::edge_link(int st, int ed, int hrank, bool reverse, int size, int othDest)
{
	_h_ed = ed;
	_h_st = st;
	_homeRank = hrank;
	_otherRank = 0;
	_reverse = reverse;
	_size = size;
	_otherDest = othDest;
	_homeDest = st + (ed-st)/2;
	_destReverse = false;
	_OthCompMode = false;
	_reArranged = false;
	_trueReverse = false;
}

edge_link::~edge_link()
{
}

bool edge_link::PointerCompare::operator()(const edge_link* l, const edge_link* r) 
{
	return *l < *r;
}

bool edge_link::operator<(const edge_link &other) const
{
	if(_OthCompMode)
	{
		if(_destReverse)
		{
			return _otherDest > other._otherDest;
		}
		else
		{
			return _otherDest < other._otherDest;
		}
	}
	else
	{
		return _homeDest < other._homeDest;
	}
	
	return 0;
}

bool edge_link::operator>(const edge_link &other) const
{
	if(_OthCompMode)
	{
		if(_destReverse)
		{
			return _otherDest < other._otherDest;
		}
		else
		{
			return _otherDest > other._otherDest;
		}
	}
	else
	{
		return _homeDest > other._homeDest;
	}
	
	return 0;
}

void edge_link::SetOtherRank(int othRank)
{
	_otherRank = othRank;
}
void edge_link::SetRearranged(bool on)
{
	_reArranged = on;
}
void edge_link::SetDestReverse(bool on)
{
	_destReverse = on;
}

int edge_link::GetRDiff()
{
	//std::cout << _size << "\t" << _homeRank << " vs " << _otherRank << "\n";
	return  abs(_homeRank - _otherRank);
}