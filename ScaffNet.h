#ifndef SCAFFNET_H
#define SCAFFNET_H
#include <map>
#include "Scaffold.h"
#include <algorithm>

class ScaffNet
{
public:
	ScaffNet(int arrSize, int minLen);
	~ScaffNet();
	void CreateScaffold(string id, int length);
	void CalcAllInScores();
	void SendUniques(vector<st_ed*> *uniques, string scaff);
	void SendRepeats(vector<st_ed*> *repeats, string scaff);
	void SendEdge(string scaff, string edge, int weight, vector<st_ed*>* links);
	void SendStats(string scaff, int single, int totalDupe);
	void AssessGenome(int genomeSize);
	int GetScaffLength(string scaff);
	int GetScaffCount();
	void SetFragChopLen(int fraglen);
	void PrintGenomeStats(char * genomeOut);
	vector<string> *ChopUpScaffold(string sid, string dna);
	Scaffold *GetScaff(string id);
	
	struct GenomeCounter {
		GenomeCounter();
		GenomeCounter(GenomeCounter *oth);
		int Size;
		int InternalDupe;
		int UniqueSeq;
		int ScaffCount;
		int ConversionSize;
	};
	
	struct SaveLink
	{
	int h_st;
	int h_ed;
	int size;
	bool reversal;
	};
	
	map<string, vector<SaveLink*>*> *ScaffLinks;
	map<string, vector<SaveLink*>*> *RepeatLinks;

private:
	void PopulateSortVec();
	GenomeCounter* InitialOutStat();
	int GetBestScaffInVector(vector<Scaffold*> *scaffs);
	int _scaffCount;
	int _currPos;
	int _fragChop;
	
	Scaffold **ScaffArray;
	std::map<string, int> ScafMap;
	
	vector<Scaffold*> _sortVec;
	vector<GenomeCounter*> _inGenStats;
	vector<GenomeCounter*> _outGenStats;
	
	vector<string> _genomeScaffs;
	vector<string> _toChopUpScaffs;
};



#endif // SCAFFNET_H
