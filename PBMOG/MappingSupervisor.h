//
//  MappingSupervisor.h
//  PBMOG
//
//  Created by malfoy on 09/01/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#ifndef __PBMOG__MappingSupervisor__
#define __PBMOG__MappingSupervisor__


#include "Utils.h"
#include "graph.h"


using namespace std;

struct path {
	string str;
	vector<uNumber> numbers;
};

class MappingSupervisor{
public:
	ofstream outFile;
	size_t offset,minSizeUnitigs,nbThreads;
	vector<string> unitigs,reads;
	unordered_map<minimizer, vector<rNumber>> min2Reads;
	unordered_set<rNumber> readused;
	size_t k,multi,H,part,k2,kgraph;
	char depthMax,nuc;
	double minJacc;
	bool mapPartAllowed,errorInKmers;
	double errorRate;
	graph G;
	mutex mutexReadReads,mutexEraseReads;
	atomic<size_t> readMapped,aligneOnPathSucess,unitigsPreMapped,bigUnitig,island,regionmapped,leftmap,rightmap,leftmapFail,rightmapFail,candidate,fail,indice;



	MappingSupervisor(const vector<string>& Iunitigs, unordered_map<minimizer, vector<rNumber>>& Iindex, size_t Ik, const vector<string>& Ireads, size_t Imulti, size_t IH, size_t Ipart, size_t Ik2, double IminJacc,graph& graphe,size_t Ikgraph){
		nuc=4;
		mapPartAllowed=false;
		errorInKmers=true;
		unitigs=Iunitigs;
		min2Reads=Iindex;
		k=Ik;
		reads=Ireads;
		multi=Imulti;
		H=IH;
		part=Ipart;
		k2=Ik2;
		minJacc=IminJacc;
		kgraph=Ikgraph;
		G=graphe;
		island=0;
		readMapped=0;
		aligneOnPathSucess=0;
		unitigsPreMapped=0;
		offset=100;
		minSizeUnitigs=100;
//		offset=8;
//		minSizeUnitigs=6;
		depthMax=5;
		bigUnitig=0;
		regionmapped=0;
		fail=candidate=leftmap=rightmap=leftmapFail=rightmapFail=0;
		nbThreads=4;
		errorRate=15;
		indice=0;
		outFile.open("zout.txt",ofstream::trunc);
	}



	void MapPart();
	bool mapUnitig(const string& unitig, const rNumber n, int& position, unordered_set<minimizer>& genomicKmers, string& read);
	void findCandidate(const string& unitig, unordered_set<minimizer>& minSet, unordered_map<rNumber,size_t>& Candidate, vector<unordered_set<minimizer>>& read2Min);
	void MapAll();
	bool isCandidateCorrect(const string& unitig, const string& read, unordered_set<minimizer>& genomicKmers,int& position, unordered_set<minimizer>& setmin);
	bool isCandidateCorrect2(const string& unitig, const string& read, unordered_set<minimizer>& genomicKmers,int& position, unordered_set<minimizer>& setmin);
	vector<path> listPath(size_t lengthRequired, uNumber ind, set<uNumber> usedUnitigs);

	bool alignOnPath(const path& path, const string& read, size_t position,set<uNumber>& usedUnitigs);

	vector<path> listPathSons(size_t lengthRequired, const string& substr,int depth);
	vector<path> listPathFathers(size_t lengthRequired, const string& substr, int depth);

	bool alignOnPathSons(const path& path, const string& read, size_t position,vector<uNumber>& n);
	bool alignOnPathSonsErrors(const path& path, const string& read, size_t position,vector<uNumber>& n);
	bool alignOnPathFathers(const path& path, const string& read, size_t position,vector<uNumber>& n);

	bool alignOnPathSons2(const path& path, const string& read, size_t position,vector<uNumber>&);
	bool alignOnPathFathers2(const path& path, const string& read, size_t position,vector<uNumber>&);

	string getPathEnd(vector<uNumber>& numbers);
	string getPathBegin(vector<uNumber>& numbers);

	void MapFromUnitigs(const string& unitig);
	void MapFromUnitigsErrors(const string& unitig);

	bool isCandidateCorrectMap(const string& unitig, const string& read, unordered_multimap<string,string>& genomicKmers,int& position, unordered_set<minimizer>& setMin);

	bool alignOnPathsSonsErrors(const vector<path>& path, const string& read, size_t position,vector<uNumber>& numbers);
	bool alignOnPathsSons(const vector<path>& Paths, const string& read, size_t position,vector<uNumber>& numbers);

};



#endif /* defined(__PBMOG__MappingSupervisor__) */
