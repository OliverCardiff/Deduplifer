#include<string>
#include<vector>
#include "st_ed.h"

using namespace std;

class Scaffold
{
	
public:
	Scaffold(string ID, int length);
	~Scaffold();
	bool InGenome;
	bool InCore;
	void AddEdge(Scaffold* scaff, int Edge, vector<st_ed*>* links);
	void SetCoverages(int singular, int total, int one);
	int GetUniqueSeq();
	int GetSize();
	int GetConvertSize();
	void CalcGenomeInScore();
	void PropagateTransferConsequences();
	void SetUniques(vector<st_ed*> *uni);
	void SetRepeats(vector<st_ed*> *repeats);
	void SetStats(int single, int tots);
	int GetTotalDupes(bool InGenome);
	void ToggleNetworkMode(bool on);
	void KillChances();
	void FinaliseUniques();
	double GetGenomeInScore();
	vector<string> *ChopUp(string DNA, int flen);
	string GetID();
	bool HasUniques();
	
	static int MinLen;
	static int FragLen;

	bool operator<(const Scaffold &other) const;
	bool operator>(const Scaffold &other) const;
	
	struct PointerCompare {
      bool operator()(const Scaffold* l, const Scaffold* r);
	};
	
private:
	string _mainID;
	vector<Scaffold*> *_edges;
	vector<int> *_edgeWeights;
	vector<st_ed*> *_uniques;
	vector<st_ed*> *_repeats;
	vector<vector<st_ed*>*> _edgeLinks;
	void UpdateGScore(int weight);
	void CalcConvertSize();
	void CalcConvertSizeA();
	void CalcConvertSizeB();
	vector<string> *ChopUpA(string DNA, int flen);
	vector<string> *ChopUpB(string DNA, int flen);

	int _length;
	int _singleCov;
	int _totalDupe;
	int _oneScore;
	int _edgeCount;
	int _uniqueSeq;
	int _convertSize;
	bool _hasUniques;
	bool _hasRepeats;
	
	double _genomeInScore;
	
	bool _networkMode;
};


