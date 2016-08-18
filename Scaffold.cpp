#include "Scaffold.h"
#include <iostream>
#include "Utility.h"
#include <algorithm>

using namespace std;

Scaffold::Scaffold(string ID, int length)
{
	_mainID = ID;
	_edges = new vector<Scaffold*>();
	_edgeWeights = new vector<int>();
	InGenome = false;
	_length = length;
	_singleCov = 0;
	_totalDupe = 0;
	_oneScore = 0;
	_genomeInScore = 0;
	_edgeCount = 0;
	_uniqueSeq = 0;
	_convertSize = _length;
	_networkMode = false;
	_hasUniques = false;
	_hasRepeats = false;
	InCore = false;
}

int Scaffold::MinLen = 300;
int Scaffold::FragLen = 300;

Scaffold::~Scaffold()
{
	delete _edges;
	delete _edgeWeights;
}

void Scaffold::AddEdge(Scaffold* scaff, int edge, vector<st_ed*>* links)
{
	_edges->push_back(scaff);
	_edgeWeights->push_back(edge);
	_edgeCount++;
	_edgeLinks.push_back(links);
}

void Scaffold::SetCoverages(int singular, int total, int one)
{
	_singleCov = singular;
	_totalDupe = total;
	_oneScore = one;
}
	
int Scaffold::GetUniqueSeq()
{
	_uniqueSeq = _length - _singleCov;
	return _convertSize;
}

void Scaffold::CalcGenomeInScore()
{
	_genomeInScore = 0;
	int dupes = 0;
	for(int i = 0; i < _edgeCount; i++)
	{
		if((*_edges)[i]->InGenome)
		{
			dupes += (*_edgeWeights)[i];
		}
	}
	if(dupes != 0)
	{
		_genomeInScore = (double)(GetUniqueSeq()) / (double)dupes;
	}
	else
	{
		_genomeInScore = GetUniqueSeq();
	}
}

void Scaffold::PropagateTransferConsequences()
{
	for(int i = 0; i < _edgeCount; i++)
	{
		if(!(*_edges)[i]->InGenome)
		{
			(*_edges)[i]->UpdateGScore((*_edgeWeights)[i]);
		}
	}
}

void Scaffold::UpdateGScore(int weight)
{
	if(GetUniqueSeq() != 0)
	{
		double dupeScore = (double)_genomeInScore / (double)GetUniqueSeq();
		
		dupeScore += (double)weight;
		
		_genomeInScore = (double)(GetUniqueSeq())/dupeScore;
	}
	else
	{
		_genomeInScore = 0;
	}
}

int Scaffold::GetConvertSize()
{
	if(_hasUniques)
	{
		CalcConvertSize();
	}
	else
	{
		_convertSize = _length;
	}

	return _convertSize;
}

void Scaffold::CalcConvertSize()
{
	double ratio = _length;
	if(_singleCov != 0)
	{
		ratio = (double)_length / (double)_singleCov;
	}
	
	if(ratio < 2)
	{
		CalcConvertSizeB();
	}
	else
	{
		CalcConvertSizeA();
	}
}

void Scaffold::CalcConvertSizeA()
{
	if(_hasUniques)
	{
		vector<st_ed*> postDup;
		_convertSize = 0;
		int cumu = 0;
		
		int s = 0;
		int e = 0;
		int oldEnd = 0;
		
		bool sSet = false;
		bool eSet = false;

		st_ed *st;
		st_ed *ed;
		st_ed fStart;
		st_ed fEnd;
		//cout << "Calc convert \t" + _mainID + "\t uniques \t" <<  _uniques->size() << " \n";
		for(vector<st_ed*>::size_type i = 0; i < _uniques->size(); i++)
		{
			if((*_uniques)[i]->_ids == 's')
			{
				st = (*_uniques)[i];
				s = st->_position;
				if(i == 0)
				{
					if(s - oldEnd < MinLen)
					{
						fStart._ids = 's';
						fStart._position = 0;
						postDup.push_back(&fStart);
					}
					else
					{
						postDup.push_back(st);
					}
				}
				else if (s - oldEnd > MinLen)
				{
					postDup.push_back(ed);
					postDup.push_back(st);
				}
			}
			if((*_uniques)[i]->_ids == 'e')
			{
				ed = (*_uniques)[i];
				e = ed->_position;
				oldEnd = e;
				
				if(i == _uniques->size() - 1)
				{
					if(_length - e < MinLen)
					{
						fEnd._ids = 'e';
						fEnd._position = _length;
						postDup.push_back(&fEnd);
					}
					else
					{
						postDup.push_back(ed);
					}
				}
			}
		}
		//cout << "Calc convert \t" + _mainID + "\t postDup \t" <<  postDup.size() << " \n";
		for(vector<st_ed*>::size_type i = 0; i < postDup.size(); i++)
		{
			if(postDup[i]->_ids == 's')
			{
				s = postDup[i]->_position;
				sSet = true;
				
			}
			if(postDup[i]->_ids == 'e')
			{
				e = postDup[i]->_position;
				eSet = true;
				
			}
			
			if(sSet && eSet)
			{
				if(e-s > FragLen)
				{
					cumu += e-s;
				}

				sSet = false;
				eSet = false;
			}
		}
		_convertSize = cumu;
		_uniqueSeq = _length - _singleCov;
	}
	else
	{
		_uniqueSeq = _length - _singleCov;
		_convertSize = _length - _singleCov;
	}
}

void Scaffold::CalcConvertSizeB()
{
	if(_hasUniques)
	{
		vector<st_ed*> postDup;
		vector<st_ed*> finalDup;
		_convertSize = 0;
		int cumu = 0;
		
		int s = 0;
		int e = 0;
		int oldEnd = 0;
		
		bool sSet = false;
		bool eSet = false;

		st_ed *st;
		st_ed *ed;
		st_ed fStart;
		st_ed fEnd;
		
		for(vector<st_ed*>::size_type i = 0; i < _uniques->size(); i++)
		{
			if((*_uniques)[i]->_ids == 's')
			{
				s = (*_uniques)[i]->_position;
				st = (*_uniques)[i];
				sSet = true;
				
			}
			if((*_uniques)[i]->_ids == 'e')
			{
				e = (*_uniques)[i]->_position;
				ed = (*_uniques)[i];
				eSet = true;
				
			}
			
			if(sSet && eSet)
			{
				if(e-s > FragLen)
				{
					postDup.push_back(st);
					postDup.push_back(ed);
				}

				sSet = false;
				eSet = false;
			}
		}
		
		for(vector<st_ed*>::size_type i = 0; i < postDup.size(); i++)
		{
			if(postDup[i]->_ids == 's')
			{
				st = postDup[i];
				s = st->_position;
				if(i == 0)
				{
					if(s - oldEnd < MinLen)
					{
						fStart._ids = 's';
						fStart._position = 0;
						finalDup.push_back(&fStart);
					}
					else
					{
						finalDup.push_back(st);
					}
				}
				else if (s - oldEnd > MinLen)
				{
					finalDup.push_back(ed);
					finalDup.push_back(st);
				}
			}
			if(postDup[i]->_ids == 'e')
			{
				ed = postDup[i];
				e = ed->_position;
				oldEnd = e;
				
				if(i == _uniques->size() - 1)
				{
					if(_length - e < MinLen)
					{
						fEnd._ids = 'e';
						fEnd._position = _length;
						finalDup.push_back(&fEnd);
					}
					else
					{
						finalDup.push_back(ed);
					}
				}
			}
		}
		//cout << "Calc convert \t" + _mainID + "\t postDup \t" <<  postDup.size() << " \n";
		for(vector<st_ed*>::size_type i = 0; i < finalDup.size(); i++)
		{
			if(finalDup[i]->_ids == 's')
			{
				s = finalDup[i]->_position;
				sSet = true;
				
			}
			if(finalDup[i]->_ids == 'e')
			{
				e = finalDup[i]->_position;
				eSet = true;
				
			}
			
			if(sSet && eSet)
			{
				cumu += e-s;
				
				sSet = false;
				eSet = false;
			}
		}
		_convertSize = cumu;
		_uniqueSeq = _length - _singleCov;
	}
	else
	{
		_uniqueSeq = _length - _singleCov;
		_convertSize = _length - _singleCov;
	}
}

void Scaffold::SetUniques(vector<st_ed*> *uni)
{
	_uniques = uni;
	_hasUniques = true;
	CalcConvertSize();
}

void Scaffold::FinaliseUniques()
{
	vector<st_ed*> finalV;
	
	if(_edges->size() > 0)
	{
		bool cancel = true;
		for(vector<vector<st_ed*>*>::size_type i = 0; i < _edges->size(); i++)
		{
			if((*_edges)[i]->InGenome)
			{
				Utility::MergeVectors(&finalV, _edgeLinks[i]);
				cancel = false;
			}
		}
		
		if(!cancel)
		{
			vector<st_ed*> *flats;
			if(!cancel)
			{
				flats = Utility::FlattenMatches(&finalV);
				sort(flats->begin(), flats->end(), st_ed::PointerCompare());
			}
			else
			{
				flats = new vector<st_ed*>();
			}
			vector<st_ed*> *fin = flats;
			vector<st_ed*> *flat2;
			bool madeF2 = false;
			
			if(_hasRepeats)
			{
				Utility::MergeRepeatsIn(flats, _repeats, _length, MinLen);
				sort(flats->begin(), flats->end(), st_ed::PointerCompare());
				flat2 = Utility::FlattenMatches(flats);
				fin = flat2;
				madeF2 = true;
			}
			
			if(_hasUniques)
				Utility::DeleteVectorContents(_uniques);
			
			_uniques = Utility::GetUniques(fin, fin->size(), _length);
			
			//if(_uniques->size() == 1)
			//{
			//	cout << "Oh no: " << _mainID << " has only 1 _uniques\n";
			//}
			
			Utility::DeleteVectorContents(flats);
			if(madeF2) Utility::DeleteVectorContents(flat2);
		}
		else
		{
			if(_hasUniques)
				Utility::DeleteVectorContents(_uniques);
				
			_hasUniques = false;
		}
	}
	else
	{
		_hasUniques = false;
	}
}

void Scaffold::SetRepeats(vector<st_ed*> *repeats)
{
	_repeats = repeats;
	//CalcConvertSize();
	_hasRepeats = true;
}

void Scaffold::SetStats(int single, int tots)
{
	_singleCov = single;
	_totalDupe = tots;
	//GetUniqueSeq();
}

void Scaffold::ToggleNetworkMode(bool on)
{
	_networkMode = on;
}

int Scaffold::GetTotalDupes(bool InGenome)
{
	int tot = 0;
	for(int i = 0; i < _edgeCount; i++)
	{
		if((*_edges)[i]->InGenome == InGenome)
		{
			tot += (*_edgeWeights)[i];
		}
	}
	
	if(tot > _singleCov)
		return _singleCov;
	else
		return tot;
}

int Scaffold::GetSize()
{
	return _length;
}

bool Scaffold::operator<(const Scaffold &other) const
{
	if(!_networkMode)
	{
		return _convertSize < other._convertSize;
	}
	else
	{
		return _genomeInScore < other._genomeInScore;
	}
}

bool Scaffold::operator>(const Scaffold &other) const
{
	if(!_networkMode)
	{
		return _convertSize > other._convertSize;
	}
	else
	{
		return _genomeInScore > other._genomeInScore;
	}
}

bool Scaffold::PointerCompare::operator()(const Scaffold* l, const Scaffold* r) 
{
	return *l < *r;
}

void Scaffold::KillChances()
{
	_genomeInScore = 0;
}

double Scaffold::GetGenomeInScore()
{
	return _genomeInScore;
}

vector<string> *Scaffold::ChopUp(string DNA, int flen)
{
	double ratio = _length;
	if(_singleCov != 0)
	{
		ratio = (double)_length / (double)_singleCov;
	}
	
	if(ratio < 2)
	{
		return ChopUpB(DNA, flen);
	}
	else
	{
		return ChopUpA(DNA, flen);
	}
}

vector<string> *Scaffold::ChopUpA(string DNA, int flen)
{
	vector<string> *vs = new vector<string>();
	
	vector<st_ed*> postDup;
	
	int s = 0;
	int e = 0;
	int oldEnd = 0;
	
	bool sSet = false;
	bool eSet = false;

	st_ed *st;
	st_ed *ed;
	st_ed *fStart = new st_ed();
	st_ed *fEnd = new st_ed();
	bool usedFS = false;
	bool usedFE = false;
	
	for(vector<st_ed*>::size_type i = 0; i < _uniques->size(); i++)
	{
		if((*_uniques)[i]->_ids == 's')
		{
			st = (*_uniques)[i];
			s = st->_position;
			if(i == 0)
			{
				if(s - oldEnd < MinLen)
				{
					fStart->_ids = 's';
					fStart->_position = 0;
					postDup.push_back(fStart);
					usedFS = true;
				}
				else
				{
					postDup.push_back(st);
				}
			}
			else if (s - oldEnd > MinLen)
			{
				postDup.push_back(ed);
				postDup.push_back(st);
			}
		}
		if((*_uniques)[i]->_ids == 'e')
		{
			ed = (*_uniques)[i];
			e = ed->_position;
			oldEnd = e;
			
			if(i == _uniques->size() - 1)
			{
				if(_length - e < MinLen)
				{
					fEnd->_ids = 'e';
					fEnd->_position = _length;
					postDup.push_back(fEnd);
					usedFE = true;
				}
				else
				{
					postDup.push_back(ed);
				}
			}
		}
	}
	
	for(vector<st_ed*>::size_type i = 0; i < postDup.size(); i++)
	{
		if(postDup[i]->_ids == 's')
		{
			s = postDup[i]->_position;
			sSet = true;
			
		}
		if(postDup[i]->_ids == 'e')
		{
			e = postDup[i]->_position;
			eSet = true;
			
		}
		
		if(sSet && eSet)
		{
			if(e-s > flen)
			{
				vs->push_back(DNA.substr(s, e-s));
			}

			sSet = false;
			eSet = false;
		}
	}
	
	if(!usedFE) delete fEnd;
	if(!usedFS) delete fStart;
	
	return vs;
	
}

vector<string> *Scaffold::ChopUpB(string DNA, int flen)
{
	vector<string> *vs = new vector<string>();
	
	vector<st_ed*> postDup;
	vector<st_ed*> finalDup;
	
	int s = 0;
	int e = 0;
	int oldEnd = 0;
	
	bool sSet = false;
	bool eSet = false;

	st_ed *st;
	st_ed *ed;
	st_ed *fStart = new st_ed();
	st_ed *fEnd = new st_ed();
	bool usedFS = false;
	bool usedFE = false;
	
	for(vector<st_ed*>::size_type i = 0; i < _uniques->size(); i++)
	{
		if((*_uniques)[i]->_ids == 's')
		{
			st = (*_uniques)[i];
			s = (*_uniques)[i]->_position;
			sSet = true;
			
		}
		if((*_uniques)[i]->_ids == 'e')
		{
			e = (*_uniques)[i]->_position;
			ed = (*_uniques)[i];
			eSet = true;
			
		}
		
		if(sSet && eSet)
		{
			if(e-s > flen)
			{
				postDup.push_back(st);
				postDup.push_back(ed);
			}

			sSet = false;
			eSet = false;
		}
	}
	
	for(vector<st_ed*>::size_type i = 0; i < postDup.size(); i++)
	{
		if(postDup[i]->_ids == 's')
		{
			st = postDup[i];
			s = st->_position;
			
			if(i == 0)
			{
				if(s - oldEnd < MinLen)
				{
					fStart->_ids = 's';
					fStart->_position = 0;
					finalDup.push_back(fStart);
					usedFS = true;
				}
				else
				{
					finalDup.push_back(st);
				}
			}
			else if (s - oldEnd > MinLen)
			{
				finalDup.push_back(ed);
				finalDup.push_back(st);
			}
		}
		if((*_uniques)[i]->_ids == 'e')
		{
			ed = postDup[i];
			e = ed->_position;
			oldEnd = e;
			
			if(i == _uniques->size() - 1)
			{
				if(_length - e < MinLen)
				{
					fEnd->_ids = 'e';
					fEnd->_position = _length;
					finalDup.push_back(fEnd);
					usedFE = true;
				}
				else
				{
					finalDup.push_back(ed);
				}
			}
		}
	}
	
	for(vector<st_ed*>::size_type i = 0; i < finalDup.size(); i++)
	{
		if(finalDup[i]->_ids == 's')
		{
			s = postDup[i]->_position;
			sSet = true;
			
		}
		if(finalDup[i]->_ids == 'e')
		{
			e = postDup[i]->_position;
			eSet = true;
			
		}
		
		if(sSet && eSet)
		{
			vs->push_back(DNA.substr(s, e-s));
			sSet = false;
			eSet = false;
		}
	}
	
	if(!usedFE) delete fEnd;
	if(!usedFS) delete fStart;
	
	return vs;
	
}

bool Scaffold::HasUniques()
{
	return _hasUniques;
}

string Scaffold::GetID()
{
	return _mainID;
}