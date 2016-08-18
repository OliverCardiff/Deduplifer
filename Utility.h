#include <vector>
#include <string>
#include <sstream>
#include "st_ed.h"
#include "edge_link.h"
#ifndef UTILITY_H
#define UTILITY_H

using namespace std;

class Utility
{
	static Utility* ms_instance;

public:
	static Utility* Instance();
	static void Release();
	static vector<string> split(const string &s, char delim);
	static void MergeVectors(vector<st_ed*> *ref, vector<st_ed*> *repeats);
	static vector<st_ed*> *FlattenMatches(vector<st_ed*> *matches);
	static void MergeRepeatsIn(vector<st_ed*> *flat, vector<st_ed*> *repeats, int scLen, int minlen);
	static vector<st_ed*>* GetUniques(vector<st_ed*> *flattened, vector<st_ed*>::size_type length, int scaffLen);
	static void DeleteVectorContents(vector<st_ed*> *todie);
	static void DeleteVectorContents(vector<edge_link*> *todie);
	static void SetRepeatIndependence(int rep);

private:
	Utility();
	~Utility();
	static vector<string> split(const string &s, char delim, vector<string> &elems);
	static int _repeatIndependence;

};

#endif // UTILITY_H
