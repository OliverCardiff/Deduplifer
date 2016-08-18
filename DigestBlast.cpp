#include "Utility.h"
#include <fstream>
#include "DigestBlast.h"
#include <string>
#include <iostream>
#include "st_ed.h"
#include <algorithm>
#include "edge_link.h"
#include <map>
#include "SaveLink.h"

using namespace std;

map<string, vector<ScaffNet::SaveLink*>*> scaffLinks;
map<string, vector<ScaffNet::SaveLink*>*> repeatLinks;
map<string, vector<st_ed*>*> scaffRepeats;

int ScanFaiForLength(string *filename, int minLen)
{
	ifstream fai;
	fai.open(filename->c_str());
	int scafCount = 0;
	string line;
	
	if(fai.is_open())
	{
		while(getline(fai, line))
		{
			if(line.length() > 5)
			{
				scafCount++;
			}
		}
	}
	else
	{
		cerr << filename->c_str() << " isn't a valid index file (run samtools faidx)\n";
		return 0;
	}
	return scafCount;
}

int ReadFaiToScaffNet(ScaffNet* net, string *filename)
{
	ifstream fai;
	fai.open(filename->c_str());
	string line;
	
	if(fai.is_open())
	{
		while(getline(fai, line))
		{
			if(line.length() > 5)
			{
				vector<string> sps = Utility::split(line, '\t');
				if(sps.size() > 1)
				{
					net->CreateScaffold(sps[0], std::atoi(sps[1].c_str()));
					//cout << "Sent " << sps[0] << " to scaff with no crash\n";
				}
			}
		}
	}
	else
	{
		cerr << filename->c_str() << " isn't a valid index file (run samtools faidx)\n";
		return 0;
	}
	
	return 1;
}

int GetTotalCoverage(vector<st_ed*> *matches)
{
	int cumu = 0;
	
	int s = 0;
	int e = 0;
	bool sSet = false;
	bool eSet = false;
	for(vector<st_ed*>::size_type i = 0; i < matches->size(); i++)
	{
		if((*matches)[i]->_ids == 's')
		{
			s = (*matches)[i]->_position;
			sSet = true;
			
		}
		if((*matches)[i]->_ids == 'e')
		{
			e = (*matches)[i]->_position;
			eSet = true;
			
		}
		
		if(sSet && eSet)
		{
			cumu += e-s;
			sSet = false;
			eSet = false;
		}
	}
	
	return cumu;
}

vector<edge_link*> *MakeEdgeLinks(vector<st_ed*> *flats)
{
	int s = 0; int e = 0;
	int ds = 0; int de = 0;
	int Hrank = 1;
	bool sSet = false;
	bool eSet = false;
	
	vector<edge_link*> *nextLinks = new vector<edge_link*>();

	for(vector<st_ed*>::size_type i = 0; i < flats->size(); i++)
	{	

		if((*flats)[i]->_ids == 's')
		{
			s = (*flats)[i]->_position;
			ds = (*flats)[i]->_dest;
			sSet = true;
		}
		if((*flats)[i]->_ids == 'e')
		{
			e = (*flats)[i]->_position;
			de = (*flats)[i]->_dest;
			eSet = true;
		}
		
		if(sSet && eSet)
		{
			bool rev;
			if(de > ds)
				rev = false;
			else
				rev = true;
			
			if(de < ds)
			{
				int temp = de;
				de = ds;
				ds = temp;
			}
			int mid = ds + ((de - ds)/2);
			nextLinks->push_back(new edge_link(s, e, Hrank, rev, e-s, mid));
			sSet = false;
			eSet = false;
			Hrank++;
		}
	}
	
	return nextLinks;
}

void SetAlignReverse(vector<edge_link*> *edgeLinks)
{
	int revCnt = 0;
	int Cnt = 0;
	
	for(vector<edge_link*>::size_type i = 0; i < edgeLinks->size(); i++)
	{
		Cnt++;
		if((*edgeLinks)[i]->_reverse)
		{
			revCnt++;
		}
	}
	
	if(revCnt > Cnt/2)
	{
		for(vector<edge_link*>::size_type i = 0; i < edgeLinks->size(); i++)
		{
			
			(*edgeLinks)[i]->_destReverse = true;
		}
	}
}

void SetEdgeLinkRanks(vector<edge_link*> *edgeLinks, bool home)
{
	int Erank = 1;
	for(vector<edge_link*>::size_type i = 0; i < edgeLinks->size(); i++)
	{	
		if(home)
		{
			(*edgeLinks)[i]->_homeRank = Erank;
			(*edgeLinks)[i]->_OthCompMode = true;
		}
		else
		{
			(*edgeLinks)[i]->_otherRank = Erank;
		}
		
		Erank++;
	}
}

int FindEdgeLinkNonMatches(vector<edge_link*> *edgeLinks)
{
	int nonMatches = 0;
	
	for(vector<edge_link*>::size_type i = 0; i < (*edgeLinks).size(); i++)
	{
		if((*edgeLinks)[i]->_homeRank != (*edgeLinks)[i]->_otherRank)
		{
			nonMatches++;
		}
	}
	
	return nonMatches;
}

int FindEdgeLinkBiggestDiff(vector<edge_link*> *edgeLinks)
{
	int diff = 0;
	int chosenIndex = 0;
	
	for(vector<edge_link*>::size_type i = 0; i < edgeLinks->size(); i++)
	{
		if((*edgeLinks)[i]->_homeRank != (*edgeLinks)[i]->_otherRank)
		{
			if(diff < (*edgeLinks)[i]->GetRDiff())
			{
				diff = (*edgeLinks)[i]->GetRDiff();
				chosenIndex = i;
			}
		}
	}
	
	return chosenIndex;
}

int GetNextAvailableSlot(bool *avails, int cnt)
{
	int nxt = 0;
	bool status = avails[0];
	
	while(status && nxt < cnt)
	{
		nxt++;
		status = avails[nxt];
	}
	return nxt;
}

int ResolveRankings(vector<edge_link*> *edgeLinks, string scaff)
{
	int rearrangements = 0;
	int nonMatches = FindEdgeLinkNonMatches(edgeLinks);
	bool avails[edgeLinks->size()];
	int avCnt = (int)edgeLinks->size();
	
	vector<edge_link*> sortLinks = *edgeLinks;
	while(nonMatches > 0)
	{
		//cout << avCnt << ": edges "<< nonMatches << ": non Matches\n";
		int chosen = FindEdgeLinkBiggestDiff(edgeLinks);

		(*edgeLinks)[chosen]->_reArranged = true;
		(*edgeLinks)[chosen]->_otherRank = (*edgeLinks)[chosen]->_homeRank;
		rearrangements++;
		
		ScaffNet::SaveLink *sv = new ScaffNet::SaveLink();
		sv->h_st = (*edgeLinks)[chosen]->_h_st;
		sv->h_ed = (*edgeLinks)[chosen]->_h_ed;
		sv->size = 50;
		sv->reversal = (*edgeLinks)[chosen]->_trueReverse;

		vector<ScaffNet::SaveLink*> *sl = scaffLinks[scaff];
		
		sl->push_back(sv);
		for(vector<edge_link*>::size_type i = 0; i < sortLinks.size(); i++)
		{
			avails[i] = false;
		}
		for(vector<edge_link*>::size_type i = 0; i < sortLinks.size(); i++)
		{
			if(sortLinks[i]->_homeRank == sortLinks[i]->_otherRank)
				avails[sortLinks[i]->_otherRank] = true;
		}
		for(vector<edge_link*>::size_type i = 0; i < sortLinks.size(); i++)
		{
			if(sortLinks[i]->_homeRank != sortLinks[i]->_otherRank)
			{
				int ind = GetNextAvailableSlot(avails, avCnt);
				avails[ind] = true;
				sortLinks[i]->_otherRank = ind + 1;
				if(ind + 1 == sortLinks[i]->_homeRank)
				{
					//cout << ind + 1 << "\n";
				}
			}
		}
		//cout << nonMatches << "\n";
		nonMatches = FindEdgeLinkNonMatches(edgeLinks);
		//cout << "here5\n";
	}
	
	//cout << sortLinks.size() << ": revs vs rearr :" << rearrangements << "\n";
	return rearrangements;
}

int CountEdgeLinkReversals(vector<edge_link*> *edgeLinks)
{
	int revs = 0;
	for(vector<edge_link*>::size_type i = 0; i < edgeLinks->size(); i++)
	{
		if((*edgeLinks)[i]->_reverse && !(*edgeLinks)[i]->_destReverse)
		{
			revs++;
			(*edgeLinks)[i]->_trueReverse = true;
		}
		else if(!(*edgeLinks)[i]->_reverse && (*edgeLinks)[i]->_destReverse)
		{
			revs++;
			(*edgeLinks)[i]->_trueReverse = true;
		}
	}
	return revs;
}

void ProcessEdgeLinks(vector<edge_link*> *edgeLinks, int cov, ofstream *edge, ofstream *links, string scaff, int length)
{
	SetAlignReverse(edgeLinks);

	if(scaffLinks.find(scaff) == scaffLinks.end())
	{
		scaffLinks[scaff]= new vector<ScaffNet::SaveLink*>(); 
	}

	int reversals = CountEdgeLinkReversals(edgeLinks);

	//Sort Based on the Query Scaff locations
	sort(edgeLinks->begin(), edgeLinks->end(), edge_link::PointerCompare());
	SetEdgeLinkRanks(edgeLinks, true);

	//Sort Based on the Reference Scaff locations
	sort(edgeLinks->begin(), edgeLinks->end(), edge_link::PointerCompare());

	SetEdgeLinkRanks(edgeLinks, false);

	int reArrs = ResolveRankings(edgeLinks, scaff);

	(*edge) << scaff << "\t" << length << "\t" << edgeLinks->size() << "\t" << cov << "\t" << reArrs << "\t" << reversals << "\n";
	for(vector<edge_link*>::size_type i = 0; i < edgeLinks->size(); i++)
	{
		(*links) << (*edgeLinks)[i]->_size << "\t" << (*edgeLinks)[i]->_reArranged << "\t" << (*edgeLinks)[i]->_trueReverse << "\n";
	}
}

int ReadBlastFile(ScaffNet *net, char *filename, int minLen, char *genomeOut)
{
	ifstream blast;
	blast.open(filename);
	ofstream edgeStats;
	ofstream linkStats;
	string baseFile(genomeOut);
	string edgeFile = baseFile + "_edgeStats.txt";
	string linkFile = baseFile + "_linkStats.txt";
	
	edgeStats.open(edgeFile.c_str());
	linkStats.open(linkFile.c_str());
	if(!(edgeStats.is_open() && linkStats.is_open()))
	{
		cerr << "Unable to open " + edgeFile + " or " + linkFile + " for editing\n";
		return 0;
	}
	edgeStats << "Scaffold\tScaff_len\tLink_count\tEdge_size\tRearranges\tReversals\n";
	linkStats << "Link_size\tRearranged\treversed\n";
	
	string line;
	string oldID = "none";
	string oldQ = "none";
	int lineCount = 0;
	
	vector<st_ed*> *HitVec;
	vector<st_ed*> *QVec;
	int RollingCoverage = 0;
	int GlobalCoverage = 0;
	HitVec = new vector<st_ed*>();
	QVec = new vector<st_ed*>();
	int edgeCnt = 0;
	int nodeCnt = 0;
	
	if(blast.is_open())
	{
		while(getline(blast, line))
		{
			lineCount++;
			if(lineCount % 100000 == 0)
			{
				//cout << "Processed " << lineCount << " lines of " << filename << "\n";
			}
			
			vector<string> sps = Utility::split(line, '\t');
			int start = atoi(sps[6].c_str());
			int end = atoi(sps[7].c_str());
			int d_start = atoi(sps[8].c_str());
			int d_end = atoi(sps[9].c_str());
			
			if(start > end)
			{
				int temp = start;
				start = end;
				end = temp;
			}
			//cout << sps[0] << ": actually got here\n";
			if(sps.size() > 1)
			{				
				if(sps[1].compare(oldID) != 0)
				{
					if(HitVec->size() > 0)
					{
						
						sort(HitVec->begin(), HitVec->end(), st_ed::PointerCompare());
						
						vector<st_ed*> *flat = Utility::FlattenMatches(HitVec);
						
						int qlen = net->GetScaffLength(oldQ);
						int rlen = net->GetScaffLength(oldID);
						
						int cov = GetTotalCoverage(flat);
						RollingCoverage += cov;
						
						if(cov > qlen/5 || cov > rlen/5 || cov > 5000)
						{
							vector<edge_link*> *edgeLinks = MakeEdgeLinks(HitVec);
							edgeCnt++;
							ProcessEdgeLinks(edgeLinks, cov, &edgeStats, &linkStats, sps[0], qlen);
							net->SendEdge(oldQ, oldID, cov, flat);
							Utility::MergeVectors(QVec, flat);
							Utility::DeleteVectorContents(edgeLinks);
						}
						//DeleteVectorContents(flat);
						Utility::DeleteVectorContents(HitVec);
						
						HitVec = new vector<st_ed*>();
					}
				}
				if(sps[0].compare(oldQ) != 0)
				{
					if(QVec->size() > 0)
					{
						nodeCnt++;
						if(nodeCnt % 500 == 0)
						{
							cout << "Processed " << nodeCnt << " scaffolds of " << filename << "\n";
						}
						
						sort(QVec->begin(), QVec->end(), st_ed::PointerCompare());
						vector<st_ed*> *flat = Utility::FlattenMatches(QVec);
						//cout << flat->size() << ": this is the size before\n";
						vector<st_ed*> *fin = flat;
						vector<st_ed*> *flat2;
						bool madeFlat2 = false;
						int scLen = net->GetScaff(oldQ)->GetSize();
						
						if(scaffRepeats.find(oldQ) != scaffRepeats.end())
						{
							net->SendRepeats(scaffRepeats[oldQ], oldQ);
							Utility::MergeRepeatsIn(flat, scaffRepeats[oldQ], scLen, minLen);
							//cout << flat->size() << ": this is the size after\n";
							sort(flat->begin(), flat->end(), st_ed::PointerCompare());
							flat2 = Utility::FlattenMatches(flat);
							//cout << flat2->size() << ": this is the size final\n";
							madeFlat2 = true;
							
							fin = flat2;
						}
						
						int cov = GetTotalCoverage(fin);
						
						vector<st_ed*> *unique = Utility::GetUniques(fin, fin->size(), net->GetScaff(oldQ)->GetSize());
						
						net->SendStats(oldQ, cov, RollingCoverage);
						net->SendUniques(unique, oldQ);
						
						GlobalCoverage += cov;
						RollingCoverage = 0;
						Utility::DeleteVectorContents(flat);
						if(madeFlat2)
							Utility::DeleteVectorContents(flat2);
						
						delete QVec;
						
						QVec = new vector<st_ed*>();
					}
				}
				oldID = sps[1];
				oldQ = sps[0];
								
				if(!blast.eof() && minLen < end - start)
				{
					st_ed *H_st = new st_ed(start, 's');
					st_ed *H_ed = new st_ed(end, 'e');
					//st_ed *Q_st = new st_ed(start, 's');
					//st_ed *Q_ed = new st_ed(end, 'e');
					H_st->_dest = d_start;
					H_ed->_dest = d_end;
					//Q_st->_dest = d_start;
					//Q_ed->_dest = d_end;
					
					HitVec->push_back(H_st);
					HitVec->push_back(H_ed);
					//QVec->push_back(Q_st);
					//QVec->push_back(Q_ed);
				}
			}
		}
	}
	else
	{
		cerr << filename << " failed to open blast results file (must be 'outfmt 6')\n";
		return 0;
	}
	cout << "Total dupe/repeat coverage: " << GlobalCoverage << " bp\n";
	cout << "Made: " << nodeCnt << " nodes and: " << edgeCnt << " edges\n";
	
	
	blast.close();
	edgeStats.close();
	linkStats.close();
	delete HitVec;
	delete QVec;
	
	net->ScaffLinks = &scaffLinks;
	
	return lineCount;
}

int ReadTBLFile(ScaffNet *net, const char *filename)
{
	ifstream tbl;
	string line;
	string oldID;
	vector<st_ed*> *hitVec = new vector<st_ed*>();
	tbl.open(filename);
	int cnt = 1;
	int scfs = 1;
	
	if(tbl.is_open())
	{
		while(getline(tbl, line))
		{
			if(cnt > 3)
			{
				vector<string> sps = Utility::split(line, ' ');
				
				if(oldID.compare(sps[4]) != 0)
				{
					sort(hitVec->begin(), hitVec->end(), st_ed::PointerCompare());
					scaffRepeats[oldID] = Utility::FlattenMatches(hitVec);
					
					repeatLinks[oldID] = new vector<ScaffNet::SaveLink*>(); 
					bool hasSt = false; bool hasEd = false;
					int stVal; int edVal;
					for(vector<st_ed*>::size_type i = 0; i < scaffRepeats[oldID]->size(); i++)
					{
						if((*scaffRepeats[oldID])[i]->_ids == 's')
						{
							hasSt = true;
							stVal = (*scaffRepeats[oldID])[i]->_position;
						}
						else if((*scaffRepeats[oldID])[i]->_ids == 'e')
						{
							hasEd = true;
							edVal = (*scaffRepeats[oldID])[i]->_position;
						}
						if(hasSt && hasEd)
						{
							if(abs(edVal - stVal) > 200)
							{
								ScaffNet::SaveLink *sv = new ScaffNet::SaveLink();
								sv->h_st = stVal;
								sv->h_ed = edVal;
								sv->size = 100;
								repeatLinks[oldID]->push_back(sv);
							}
							hasSt = false;
							hasEd = false;
						}
					}
					
					Utility::DeleteVectorContents(hitVec);
					hitVec = new vector<st_ed*>();
					scfs++;
					if(scfs % 10000 == 0)
					{
						cout << "Processed " << scfs << " scaffold repeat masks\n";
					}
				}
				int sw = atoi(sps[0].c_str());
				if(sw > 30)
				{
					hitVec->push_back(new st_ed(atoi(sps[5].c_str()), 's'));
					hitVec->push_back(new st_ed(atoi(sps[6].c_str()), 'e'));
				}
				oldID = sps[4];
			}
			cnt++;
		}
		sort(hitVec->begin(), hitVec->end(), st_ed::PointerCompare());
		scaffRepeats[oldID] = Utility::FlattenMatches(hitVec);
		
		repeatLinks[oldID] = new vector<ScaffNet::SaveLink*>(); 
		bool hasSt = false; bool hasEd = false;
		int stVal; int edVal;
		for(vector<st_ed*>::size_type i = 0; i < scaffRepeats[oldID]->size(); i++)
		{
			if((*scaffRepeats[oldID])[i]->_ids == 's')
			{
				hasSt = true;
				stVal = (*scaffRepeats[oldID])[i]->_position;
			}
			else if((*scaffRepeats[oldID])[i]->_ids == 'e')
			{
				hasEd = true;
				edVal = (*scaffRepeats[oldID])[i]->_position;
			}
			if(hasSt && hasEd)
			{
				if(abs(edVal - stVal) > 200)
				{
					ScaffNet::SaveLink *sv = new ScaffNet::SaveLink();
					sv->h_st = stVal;
					sv->h_ed = edVal;
					sv->size = 100;
					repeatLinks[oldID]->push_back(sv);
				}
				hasSt = false;
				hasEd = false;
			}
		}
		Utility::DeleteVectorContents(hitVec);
	}
	else
	{
		cout << "Couldn't open Repeat Masker Table file\n";
		return 0;
	}
	
	tbl.close();
	
	net->RepeatLinks = &repeatLinks;
	return 1;
}