#include "ScaffNet.h"
#include "SaveLink.h"
#include <map>
#include <string>

int ReadGenome(string g_in, char * g_out, ScaffNet *net);
int ChopUpUniques(ScaffNet *net, ofstream *g_out);
int OutputDNAsites(vector<ScaffNet::SaveLink*> *links, string DNA, string scId, ofstream *frags);