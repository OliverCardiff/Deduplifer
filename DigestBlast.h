#include "ScaffNet.h"
#include <cstdlib>

int ReadFaiToScaffNet(ScaffNet* net, std::string *filename);
int ScanFaiForLength(std::string *filename, int minLen);
int ReadBlastFile(ScaffNet *net, char *filename, int minLen, char *genomeOut);
int ReadTBLFile(ScaffNet *net, const char *filename);