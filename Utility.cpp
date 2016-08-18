#include "Utility.h"
#include <cstdlib> // NULL 

Utility* Utility::ms_instance = NULL;

using namespace std;

Utility::Utility()
{
}

Utility::~Utility()
{
}

Utility* Utility::Instance()
{
	if (ms_instance == NULL) {
		ms_instance = new Utility();
	}
	return ms_instance;
}

void Utility::Release()
{
	if (ms_instance) {
		delete ms_instance;
	}
	ms_instance = NULL;
}

vector<string> Utility::split(const string &s, char delim)
{
	vector<string> elems;
	split(s,delim, elems);
	return elems;
}

vector<string> Utility::split(const string &s, char delim, vector<string>& elems)
{
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim)) {
		if(item.length() > 0)
		{
			elems.push_back(item);
		}
	}
	return elems;
}

int Utility::_repeatIndependence = 1000;

void Utility::SetRepeatIndependence(int rep)
{
	_repeatIndependence = rep;
}

void Utility::MergeVectors(vector<st_ed*> *ref, vector<st_ed*> *repeats)
{
	for(vector<st_ed*>::size_type i = 0; i < repeats->size(); i++)
	{
		ref->push_back((*repeats)[i]);
	}
}

vector<st_ed*>* Utility::GetUniques(vector<st_ed*> *flattened, vector<st_ed*>::size_type length, int scaffLen)
{
	vector<st_ed*> *finals = new vector<st_ed*>();
	
	int depth = 0;
	char last = 's';
	
	if(length > 1)
	{
		for(vector<st_ed*>::size_type i = 0; i < flattened->size(); i++)
		{
			if(i == 0 && (*flattened)[i]->_position != 0)
			{
				finals->push_back(new st_ed(0,'s'));
				last = 's';
				if((*flattened)[i]->_ids == 's')
				{
					finals->push_back(new st_ed((*flattened)[i]->_position, 'e'));
					last = 'e';
				}
			}
			
			if(depth > 0 && (*flattened)[i]->_ids == 'e')
			{
				finals->push_back(new st_ed((*flattened)[i]->_position, 's'));
				last = 's';
			}
			if(depth < 1 && (*flattened)[i]->_ids == 's')
			{
				finals->push_back(new st_ed((*flattened)[i]->_position, 'e'));
				last = 'e';
			}
			
			if((*flattened)[i]->_ids == 'e')
			{
				depth--;
			}
			else
				depth++;
		}
	}
	
	if(last == 's')
	{
		finals->push_back(new st_ed(scaffLen, 'e'));
	}
	
	return finals;
}

vector<st_ed*> *Utility::FlattenMatches(vector<st_ed*> *matches)
{
	vector<st_ed*> *flattened = new vector<st_ed*>();
	
	vector<st_ed*>::size_type length = matches->size();
	
	int depth = 0;
	bool skipNextStart = false;
	
	if(length > 1)
	{
		for(vector<st_ed*>::size_type i = 0; i < matches->size(); i++)
		{
			if((*matches)[i]->_ids == 's' && depth < 1 && !((*matches)[i+1]->_ids == 'e' 
				&& (*matches)[i+1]->_position == (*matches)[i]->_position))
			{
				if(!skipNextStart)
				{
					flattened->push_back(new st_ed((*matches)[i]->_position, 's'));
				}
				else
					skipNextStart = false;
			}
			if((*matches)[i]->_ids == 'e' && depth == 1 && !((*matches)[i-1]->_ids == 's' 
				&& (*matches)[i-1]->_position == (*matches)[i]->_position))
			{
				if(i < matches->size() - 1)
				{
					if((*matches)[i+1]->_ids == 's' && abs((*matches)[i+1]->_position - (*matches)[i]->_position) < 100)
						skipNextStart = true;
				}
				if(!skipNextStart)
					flattened->push_back(new st_ed((*matches)[i]->_position, 'e'));
			}
			
			if((*matches)[i]->_ids == 's')
				depth++;
			if((*matches)[i]->_ids == 'e')
				depth--;
		}
	}
	
	return flattened;
}

void Utility::DeleteVectorContents(vector<st_ed*> *todie)
{
	if(todie->size() > 1)
	{
		for(vector<st_ed*>::size_type i = 0; i < todie->size(); i++)
		{
			delete (*todie)[i];
		}
	}
	todie->clear();
	
	delete todie;
}

void Utility::DeleteVectorContents(vector<edge_link*> *todie)
{
	if(todie->size() > 1)
	{
		for(vector<edge_link*>::size_type i = 0; i < todie->size(); i++)
		{
			delete (*todie)[i];
		}
	}
	todie->clear();
	
	delete todie;
}

void Utility::MergeRepeatsIn(vector<st_ed*> *flat, vector<st_ed*> *repeats, int scLen, int minlen)
{
	bool keeps[repeats->size()];
	for(vector<st_ed*>::size_type i = 0; i < repeats->size(); i += 2)
	{
		keeps[i] = false;
		keeps[i + 1] = false;
		if(abs((*repeats)[i]->_position - (*repeats)[i+1]->_position) > _repeatIndependence)
		{
			keeps[i] = true;
			keeps[i + 1] = true;
		}
		else if((*repeats)[i]->_position < minlen)
		{
			keeps[i] = true;
			keeps[i + 1] = true;
		}
		else if(abs(scLen - (*repeats)[i + 1]->_position) < minlen)
		{
			keeps[i] = true;
			keeps[i + 1] = true;
		}
		else
		{
			for(vector<st_ed*>::size_type j = 0; j < flat->size(); j++)
			{
				if(abs((*repeats)[i]->_position - (*flat)[j]->_position)  < minlen)
				{
					keeps[i] = true;
					keeps[i + 1] = true;
				}
				else if(abs((*repeats)[i + 1]->_position - (*flat)[j]->_position)  < minlen)
				{
					keeps[i] = true;
					keeps[i + 1] = true;
				} 
			}
		}
	}
	for(vector<st_ed*>::size_type i = 0; i < repeats->size(); i++)
	{
		if(keeps[i])
		{
			flat->push_back(new st_ed((*repeats)[i]->_position, (*repeats)[i]->_ids));
		}
	}
}