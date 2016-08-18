#include "DigestGenome.h"
#include "ScaffNet.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

map<string, string> savedScaffs;

void PrintDNALines(string DNA, int len, ofstream *ofs)
{
	string temp;
	string intermed = DNA.c_str();
	
	while(DNA.length() > len)
	{
		temp = DNA.substr(0, len);
		(*ofs) << temp << "\n";
		DNA.erase(0, len);
	}
	(*ofs) << DNA << "\n";
}

int ChopUpUniques(ScaffNet *net, ofstream *g_out)
{
	cout << "Chopping up remaining scaffolds\n";
	int totCnt = 0;
	int scaffs = 0;
	for(map<string, string>::iterator it = savedScaffs.begin(); it != savedScaffs.end(); ++it)
	{
		int cnt = 1;
		string DNA = it->second;
		string scaf = it->first;
		Scaffold *sc = net->GetScaff(scaf);
		scaffs++;
		if(sc->GetUniqueSeq() > 0)
		{
			//cout << "Chopping " << scaf << "\n";
			//This check is needed to find out if there are markings on the contig indicating chop up points!
			if(sc->HasUniques())
			{
				vector<string> *vcs = net->ChopUpScaffold(scaf, DNA);
				
				if(vcs->size() > 0)
				{
					//cout << "Chopped into: " << vcs->size() << "\n";
					
					for(vector<string>::size_type j = 0; j < vcs->size(); j++)
					{
						totCnt++;
						(*g_out) << ">" << scaf << "_chop_" << cnt << "\n";
						PrintDNALines((*vcs)[j], 100, g_out);
						cnt++;
					}
					delete vcs;
				}
			}
			else
			{
				totCnt++;
				(*g_out) << ">" << scaf << "_chop_" << cnt << "\n";
				PrintDNALines(DNA, 100, g_out);
			}
		}
	}
	
	cout << "Chopped: " << scaffs << " scaffolds into: " << totCnt << " contigs\n";
	cout << "Rescaffolding is advised!\n";
	return 1;
}

int ReadGenome(string g_in, char * g_out, ScaffNet *net)
{	
	string fileBase(g_out);
	string repsFile = fileBase + "_repeatEdges.fasta";
	string fragsFile = fileBase + "_fragments.fasta";
	string genomeOut = fileBase + ".reduced.fasta";
	ifstream genome;
	ofstream newGenome;
	ofstream fragments;
	ofstream repeats;
	repeats.open(repsFile.c_str());
	fragments.open(fragsFile.c_str());
	genome.open(g_in.c_str());
	newGenome.open(genomeOut.c_str());
	
	if(!fragments.is_open())
	{
		cerr << "Unable to open/create fragments file " + fragsFile + "\n";
		return 0;
	}
	if(!repeats.is_open())
	{
		cerr << "Unable to open/create repeat Edges file " + repsFile + "\n";
		return 0;
	}
	if(!newGenome.is_open())
	{
		cerr << "Unable to open/create new genome fasta file " + genomeOut + "\n";
		return 0;
	}
	
	bool flip = 0;
	bool chop = 0;
	string line;
	string oldCopy;
	string concat = "";
	int cnt = 0;
	if(genome.is_open())
	{
		while(getline(genome, line))
		{
			const char* cline = line.c_str();
			//cout << cline[0] << "\n";
			if(cline[0] == '>')
			{
				string copy = line.c_str();
				copy.erase(0,1);
				Scaffold *scf = net->GetScaff(copy);
				
				if(cnt != 0)
				{
					if(chop)
					{
						savedScaffs[oldCopy.c_str()] = concat.c_str();
					}
				}
				
				if(net->ScaffLinks->find(oldCopy) != net->ScaffLinks->end())
				{
					vector<ScaffNet::SaveLink*> *vec = net->ScaffLinks->at(oldCopy);
					OutputDNAsites(vec, concat, oldCopy, &fragments);
				}
				
				if(net->RepeatLinks->find(oldCopy) != net->RepeatLinks->end())
				{
					vector<ScaffNet::SaveLink*> *vec = net->RepeatLinks->at(oldCopy);
					OutputDNAsites(vec, concat, oldCopy, &repeats);
				}
				
				if(scf->InGenome && scf->InCore)
				{
					flip = true;
				}
				else
					flip = false;
					
				if(!scf->InCore && scf->InGenome)
				{
					chop = true;
				}
				else
					chop = false;
				
				concat = "";
				oldCopy = copy.c_str();
			}
			else
			{
				if(line.length() > 0)
				{
					concat += line;
				}
			}
			cnt++;
			if(cnt % 1000000 == 0)
			{
				cout << "read: " << cnt << " lines\n";
			}
			if(flip)
			{
				newGenome << line << "\n";
			}
		}
	}
	else
	{
		cout << "Failed to Open Genome: " + g_in + "\n";
		return 0;
	}
	
	if(!flip)
	{
		savedScaffs[oldCopy.c_str()] = concat.c_str();
	}
	
	cout << "read : " << cnt << "\n";
	genome.close();
	fragments.close();
	ChopUpUniques(net, &newGenome);
	newGenome.close();
	return 1;
}

int OutputDNAsites(vector<ScaffNet::SaveLink*> *links, string DNA, string scId, ofstream *frags)
{
	int cnt = 1;
	
	vector<ScaffNet::SaveLink*> nonOverlaps;
	
	for(vector<ScaffNet::SaveLink*>::size_type i = 0; i < links->size(); i++)
	{
		int st = (*links)[i]->h_st;
		int ed = (*links)[i]->h_ed;
		if(ed - st > 300)
		{
			if(nonOverlaps.size() > 0)
			{
				bool add = true;
				for(vector<ScaffNet::SaveLink*>::size_type j = 0; j < nonOverlaps.size(); j++)
				{	
					int st2 = (*links)[i]->h_st;
					int ed2 = (*links)[i]->h_ed;
				
					if(abs(st2 - st) < 100 || abs(ed2 - ed) < 100)
					{
						add = false;
					}
				}
				if(add)
				{
					nonOverlaps.push_back((*links)[i]);
				}
			}
			else
			{
				nonOverlaps.push_back((*links)[i]);
			}
		}
	}
	for(vector<ScaffNet::SaveLink*>::size_type i = 0; i < nonOverlaps.size(); i++)
	{
		int len = DNA.length(); 
		int st = nonOverlaps[i]->h_st;
		int ed = nonOverlaps[i]->h_ed;
		int sz = ed - st;
		
		int s1 = st - nonOverlaps[i]->size;
		int s2 = st + nonOverlaps[i]->size;
		int e1 = ed - nonOverlaps[i]->size;
		int e2 = ed + nonOverlaps[i]->size;
		
		if(s1 < 0) s1 = 0;
		if(s2 >= len) s2 = len - 1;
		if(e1 < 0) e1 = 0;
		if(e2 >= len) e2 = len - 1;
		
		if(s1 < s2)
		{
			string cut = DNA.substr(s1, s2-s1);
			(*frags) << ">" << scId << "_" << cnt << "_" << nonOverlaps[i]->reversal << "\n";
			(*frags) << cut << "\n";
			cnt++;
		}
		if(e1 < e2)
		{
			string cut = DNA.substr(e1, e2-e1);
			(*frags) << ">" << scId << "_" << cnt << "_" << nonOverlaps[i]->reversal << "\n";
			(*frags) << cut << "\n";
			cnt++;
		}
	}
	
	return 1;
}