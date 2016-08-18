#ifndef EDGE_LINK_H
#define EDGE_LINK_H


class edge_link
{
public:
	edge_link(int st, int ed, int hrank, bool reverse, int size, int othDest);
	~edge_link();
	
	bool operator<(const edge_link &other) const;
	bool operator>(const edge_link &other) const;
	
	void SetOtherRank(int othRank);
	void SetRearranged(bool on);
	void SetDestReverse(bool on);
	int GetRDiff();
	
	int _h_st;
	int _h_ed;
	int _homeRank;
	int _homeDest;
	int _otherRank;
	int _otherDest;
	int _size;
	bool _reverse;
	bool _destReverse;
	bool _reArranged;
	bool _OthCompMode;
	bool _trueReverse;
	
	struct PointerCompare {
      bool operator()(const edge_link* l, const edge_link* r);
	};
};

#endif // EDGE_LINK_H
