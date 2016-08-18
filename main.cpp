#include "DigestBlast.h"
#include <stdio.h>
#include <getopt.h>
#include <iostream>
#include "SaveLink.h"
#include "DigestGenome.h"
#include "Utility.h"

using namespace std;

bool RunTest(bool gS, bool bS, bool mS, bool oS, bool lS, bool sS, bool rS, bool fS)
{
	bool pass = true;
	if(!rS)
	{
		cout << "-r <integer> Repeat Independence: (optional) Argument absent!\n";
		pass = true;
	}
	if(!rS)
	{
		cout << "-f <integer> Fragment Size: (optional) Argument absent!\n";
		pass = true;
	}
	if(!gS)
	{
		cout << "-g <filename> Input Genome: Argument absent!\n";
		pass = false;
	}
	if(!bS)
	{
		cout << "-b <filename> blast.oufmt6: Argument absent!\n";
		pass = false;
	}
	if(!mS)
	{
		cout << "-m <filename> Repeat Masker Annotation (fasta.out): Argument absent!\n";
		pass = false;
	}
	if(!oS)
	{
		cout << "-o <base_filename> Output prefix for all files\n";
		pass = false;
	}
	if(!lS)
	{
		cout << "-l <int> Min Duplication Length: Argument absent!\n";
		pass = false;
	}
	if(!sS)
	{
		cout << "-s <int> Genome size: Argument absent!\n";
		pass = false;
	}
}

int main(int argc, char **argv)
{
	char *genomeFile;
	char *genomeOut;
	char *repeatMask;
	string *genomeIndex;
	char *blastFile;
	int minLen = 100;
	int GenomeSize = 850000000;
	int repeatInd = 5000;
	int fragSize = 500;
	ScaffNet *mainNet;
	int c;
	
	bool gS = false; bool bS = false; bool mS = false; bool rS = false;
	bool oS = false; bool lS = false; bool sS = false; bool fS = false;
	
	while((c = getopt (argc, argv, "g:b:l:p:s:o:m:r:f:")) != -1)
	{
		switch(c)
		{
			case 'g':
				gS = true;
				genomeFile = optarg;
				genomeIndex = new string(genomeFile);
				*genomeIndex = *genomeIndex + ".fai";
				break;
			case 'b':
				bS = true;
				blastFile = optarg;
				break;
			case 'f':
				fS = true;
				fragSize = atoi(optarg);
				break;
			case 'm':
				mS = true;
				repeatMask = optarg;
				break;
			case 'o':
				oS = true;
				genomeOut = optarg;
				break;
			case 'l':
				lS = true;
				minLen = atoi(optarg);
				break;
			case 'r':
				rS = true;
				repeatInd = atoi(optarg);
				break;
			case 's':
				sS = true;
				GenomeSize = atoi(optarg);
				break;
			case '?':
				if(optopt == 'g' || optopt == 'b' || optopt == 'c' || optopt == 's' || optopt == 'm' || optopt == 'o')
				{
					fprintf(stderr, "Option -%c requires an argument.\n", optopt);
				}
				else if(isprint (optopt))
				{
					fprintf(stderr, "Unknown Option '-%c'.\n", optopt);
				}
				else
				{
					fprintf(stderr, "Unknown Option '\\x%x'.\n", optopt);
				}
				return 1;
			default:
				abort();
		}
	}
	
	if(!RunTest(gS, bS, mS, oS, lS, sS, rS, fS))
	{
		return 1;
	}
	Utility::SetRepeatIndependence(repeatInd);
	int len = ScanFaiForLength(genomeIndex, minLen);
	mainNet = new ScaffNet(len, minLen);
	if(!fS)
	{
		fragSize = minLen;
	}
	mainNet->SetFragChopLen(fragSize);
	cout << "Initialising DupeNetwork, size: " << len <<"...\n";
	if(!ReadFaiToScaffNet(mainNet, genomeIndex))
		return 1;
	cout << "Reading Repeat Mask\n";
	if(!ReadTBLFile(mainNet, repeatMask))
		return 1;
	cout << "Running allele rearrangement analysis\n";
	if(!ReadBlastFile(mainNet, blastFile, minLen, genomeOut))
		return 1;
	cout << "Network Built. Model Genome Reconstruction at: " << GenomeSize << " bases\n";
	mainNet->AssessGenome(GenomeSize);
	cout << "Assessment complete, printing construction data..\n";
	mainNet->PrintGenomeStats(genomeOut);
	cout << "Reading Genome\n";
	if(!ReadGenome(genomeFile, genomeOut, mainNet))
		return 1;
	return 0;
}
