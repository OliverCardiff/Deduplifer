#include "ScaffNet.h"
#include <fstream>
#include <iostream>

using namespace std;

ScaffNet::ScaffNet(int arrSize, int minLen)
{
	_scaffCount = arrSize;
	ScaffArray = new Scaffold*[arrSize];
	_currPos = 0;
	Scaffold::MinLen = minLen;
}

ScaffNet::~ScaffNet()
{
	_sortVec.clear();
	for(int i = 0; i < _scaffCount; i++)
	{
		delete ScaffArray[i];
	}
}

ScaffNet::GenomeCounter::GenomeCounter()
{
	Size = 0;
	InternalDupe = 0;
	UniqueSeq = 0;
	ScaffCount = 0;
	ConversionSize = 0;
}

ScaffNet::GenomeCounter::GenomeCounter(GenomeCounter *oth)
{
	Size = oth->Size;
	InternalDupe = oth->InternalDupe;
	UniqueSeq = oth->UniqueSeq;
	ScaffCount =  oth->ScaffCount;
	ConversionSize =  oth->ConversionSize;
}

void ScaffNet::SetFragChopLen(int flen)
{
	_fragChop = flen;
	Scaffold::FragLen = flen;
}

void ScaffNet::CreateScaffold(string sid, int length)
{
	ScaffArray[_currPos] = new Scaffold(sid, length);
	
	ScafMap.insert(pair<string, int>(sid, _currPos));
	_currPos++;
}

void ScaffNet::CalcAllInScores()
{
	for(int i = 0; i < _scaffCount; i++)
	{
		if(!ScaffArray[i]->InGenome)
		{
			ScaffArray[i]->CalcGenomeInScore();
		}
	}
}

int ScaffNet::GetScaffCount()
{
	return _scaffCount;
}

void ScaffNet::SendUniques(vector<st_ed*> *uni,  string scaff)
{
	int ind = ScafMap[scaff];
	//cout << scaff << " = sent string\n";
	//cout << ind << " = found index\n";
	ScaffArray[ind]->SetUniques(uni);
}

void ScaffNet::SendRepeats(vector<st_ed*> *repeats, string scaff)
{
	int ind = ScafMap[scaff];
	//cout << scaff << " = sent string\n";
	//cout << ind << " = found index\n";
	ScaffArray[ind]->SetRepeats(repeats);
}

void ScaffNet::SendEdge(string scaff, string edge, int weight, vector<st_ed*>* links)
{
	int refind = ScafMap[scaff];
	int edind = ScafMap[edge];
	
	ScaffArray[refind]->AddEdge(ScaffArray[edind], weight, links);
}

void ScaffNet::SendStats(string scaff, int single, int totalDupe)
{
	int refind = ScafMap[scaff];
	
	ScaffArray[refind]->SetStats(single, totalDupe);
}

vector<string> *ScaffNet::ChopUpScaffold(string scaff, string dna)
{
	int refind = ScafMap[scaff];
	
	return ScaffArray[refind]->ChopUp(dna, _fragChop);
}

void ScaffNet::PopulateSortVec()
{
	for(int i = 0; i < _scaffCount; i++)
	{
		_sortVec.push_back(ScaffArray[i]);
	}
}

ScaffNet::GenomeCounter* ScaffNet::InitialOutStat()
{
	GenomeCounter * gen = new GenomeCounter();
	gen->ScaffCount = _scaffCount;
	gen->InternalDupe = 0;
	gen->ConversionSize = 0;
	gen->Size = 0;
	int cnt = 0;
	
	for(int i = 0; i < _scaffCount; i++)
	{
		gen->ConversionSize += ScaffArray[i]->GetUniqueSeq();
		gen->UniqueSeq += ScaffArray[i]->GetUniqueSeq();
		gen->Size += ScaffArray[i]->GetSize();
		gen->InternalDupe += ScaffArray[i]->GetTotalDupes(false);
		
		if(gen->InternalDupe < 0 && cnt < 10)
		{
			cnt++;
			cout << ScaffArray[i]->GetUniqueSeq() << " > " << ScaffArray[i]->GetSize() << "\n";
			cout << ScaffArray[i]->GetTotalDupes(false) << ": dupes\n";
			cout << ScaffArray[i]->GetID() << ": ID\n";
			cout << gen->InternalDupe << ": total \n";
		}
	}
	
	gen->InternalDupe = gen->InternalDupe / 2;
	
	return gen;
}

int ScaffNet::GetScaffLength(string scaff)
{
	int refind = ScafMap[scaff];
	return ScaffArray[refind]->GetSize();
}

int ScaffNet::GetBestScaffInVector(vector<Scaffold*> *scaffs)
{
	int ind = 0;
	double score = 0;
	
	for(vector<Scaffold*>::size_type i = 0; i < (*scaffs).size(); i++)
	{
		if((*scaffs)[i]->GetGenomeInScore() > score)
		{
			ind = i;
			score = (*scaffs)[i]->GetGenomeInScore();
		}
	}
	
	return ind;
}

void ScaffNet::AssessGenome(int genomeSize)
{
	PopulateSortVec();
	
	sort(_sortVec.begin(), _sortVec.end(), Scaffold::PointerCompare());
	
	int stage1 = genomeSize / 25;
	int currSize = 0;
	int j = _scaffCount - 1;
	int added = 0;
	GenomeCounter *prevIn = new GenomeCounter();
	GenomeCounter *prevOut = InitialOutStat();
	
	int totalScaffs = prevOut->ScaffCount;
	//cout << "Entering First Loop\n";
	while(currSize < stage1 && j > 1)
	{
		_sortVec[j]->InGenome = true;
		_sortVec[j]->InCore = true;
		int len = _sortVec[j]->GetSize();
		
		GenomeCounter *inStat = new GenomeCounter(prevIn);
		GenomeCounter *outStat = new GenomeCounter(prevOut);
		
		outStat->ConversionSize += len - _sortVec[j]->GetUniqueSeq();
		outStat->InternalDupe -= _sortVec[j]->GetTotalDupes(false)/2;
		outStat->ScaffCount--;
		outStat->Size -= len;
		outStat->UniqueSeq -= _sortVec[j]->GetUniqueSeq();
		
		inStat->InternalDupe += _sortVec[j]->GetTotalDupes(true);
		inStat->ScaffCount++;
		inStat->Size += len;
		inStat->UniqueSeq += _sortVec[j]->GetUniqueSeq();
		
		_outGenStats.push_back(outStat);
		_inGenStats.push_back(inStat);
		
		prevIn = inStat;
		prevOut = outStat;
		
		added++;
		currSize += len;
		j--;
	}
	//cout << "Exiting First Loop\n";
	cout << "Added: " << added << ", scaffolds, size: " << currSize << ", stage 1 target: " << stage1 << ", consize: " << prevOut->ConversionSize << "\n";
	CalcAllInScores();
	
	vector<Scaffold*> remains;
	
	for(int i = 0; i < _scaffCount; i++)
	{
		if(!_sortVec[i]->InGenome)
		{
			remains.push_back(_sortVec[i]);
		}
		_sortVec[i]->ToggleNetworkMode(true);
	}
	
	//int remSize = (int)remains.size();
	cout << "Performing duplication network node exchange..\n";
	while(prevIn->Size < genomeSize && added < totalScaffs - 1)
	{
		//sort(remains.begin(), remains.end(), Scaffold::PointerCompare());
		int bestInd = GetBestScaffInVector(&remains);
		
		remains[bestInd]->InGenome = true;
		remains[bestInd]->PropagateTransferConsequences();
		remains[bestInd]->KillChances();
		remains[bestInd]->FinaliseUniques();
		int len = remains[bestInd]->GetConvertSize();
		//int len = remains[bestInd]->GetSize();
		
		GenomeCounter *inStat = new GenomeCounter(prevIn);
		GenomeCounter *outStat = new GenomeCounter(prevOut);
		
		outStat->ConversionSize += len - remains[bestInd]->GetUniqueSeq();
		outStat->InternalDupe -= remains[bestInd]->GetTotalDupes(false)/2;
		outStat->ScaffCount--;
		outStat->Size -= len;
		outStat->UniqueSeq -= remains[bestInd]->GetUniqueSeq();
		
		inStat->InternalDupe += remains[bestInd]->GetTotalDupes(true);
		inStat->ScaffCount++;
		inStat->Size += len;
		inStat->UniqueSeq += remains[bestInd]->GetUniqueSeq();
		
		_outGenStats.push_back(outStat);
		_inGenStats.push_back(inStat);
		
		prevIn = inStat;
		prevOut = outStat;
		added++;
		
		if(added % 1000 == 0)
		{
			cout << "Added: " << added << " so far..\t " << "size:\t" << inStat->Size << " \n";
		} 
	}
	cout << "added: " << added << ", totScaff: " << totalScaffs << ", genom:, " << prevIn->Size << " consize: " << prevOut->ConversionSize << "\n";
	cout << "Exiting Second Loop\n";
}

Scaffold *ScaffNet::GetScaff(string id)
{
	int ind = ScafMap[id];
	return ScaffArray[ind];
}

void ScaffNet::PrintGenomeStats(char * genomeOut)
{
	string fileBase(genomeOut);
	string statFile = fileBase + "_genomeStats.txt";
	ofstream gstat;
	gstat.open(statFile.c_str());
	
	gstat << "Con_Size\tOut_Dupe\tOut_Count\tOut_Size\tOut_Unique\tIn_Dupe\tIn_Count\tIn_Size\tIn_Unique\n";
	
	for(vector<GenomeCounter*>::size_type i = 0; i < _outGenStats.size(); i++)
	{
		GenomeCounter *g1 = _outGenStats[i];
		GenomeCounter *g2 = _inGenStats[i];
		
		gstat << g1->ConversionSize << "\t" << g1->InternalDupe << "\t" << g1->ScaffCount << "\t" << g1->Size << "\t" << g1->UniqueSeq << "\t";
		gstat << g2->InternalDupe << "\t" << g2->ScaffCount << "\t" << g2->Size << "\t" << g2->UniqueSeq << "\n";
	}
	
	gstat.close();
}
