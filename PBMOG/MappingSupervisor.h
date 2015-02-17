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
	uNumber lastUnitig;
};

class MappingSupervisor{
public:
	ofstream outFile;
	size_t offset,minSizeUnitigs;
	vector<string> unitigs,reads;
	unordered_map<minimizer, unordered_set<rNumber>> min2Reads;
	size_t k,multi,H,part,k2,kgraph;
	char depthMax;
	double minJacc;
	graph G;
	mutex myMutex,mutexEraseReads;
	atomic<size_t> readMapped,aligneOnPathSucess,unitigsPreMapped,bigUnitig,island;



	MappingSupervisor(const vector<string>& Iunitigs, unordered_map<minimizer, unordered_set<rNumber>>& Iindex, size_t Ik, const vector<string>& Ireads, size_t Imulti, size_t IH, size_t Ipart, size_t Ik2, double IminJacc,graph& graphe,size_t Ikgraph){
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
		minSizeUnitigs=30;
		depthMax=5;
		bigUnitig=0;
		outFile.open("zout.txt",ofstream::trunc);
	}



	void MapPart(size_t L, size_t R);
	void findCandidate(const string& unitig, unordered_set<minimizer>& min, unordered_map<rNumber,size_t>& count, unordered_map<rNumber,unordered_set<minimizer>>& read2min);
	void MapAll();
	int isCandidateCorrect(const string& unitig, rNumber readNumber, unordered_map<rNumber,unordered_set<minimizer>>& read2min, unordered_set<minimizer>& genomicKmers);
	
	vector<path> listPath(size_t lengthRequired, uNumber ind, set<uNumber> usedUnitigs);

		vector<path> listPathSons(size_t lengthRequired, uNumber ind, unordered_set<uNumber>& usedUnitigs);
		vector<path> listPathFathers(size_t lengthRequired, uNumber ind, unordered_set<uNumber>& usedUnitigs);
	bool alignOnPath(const path& path, const string& read, size_t position,set<uNumber>& usedUnitigs);

	vector<path> listPathSons(size_t lengthRequired, uNumber ind,char);
	vector<path> listPathFathers(size_t lengthRequired, uNumber ind, char depth);

	bool alignOnPathSons(const path& path, const string& read, size_t position);
		bool alignOnPathFathers(const path& path, const string& read, size_t position);


};



#endif /* defined(__PBMOG__MappingSupervisor__) */
