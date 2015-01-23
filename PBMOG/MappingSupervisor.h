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


using namespace std;

struct path {
	string str;
	uNumber lastUnitig;
};

class MappingSupervisor{
public:
	size_t offset,minSizeUnitigs;
	vector<string> unitigs,reads;
	unordered_map<minimizer, unordered_set<rNumber>> min2Reads;
	size_t k,multi,H,part,k2,kgraph;
	double minJacc;
	unordered_map<string, vector<uNumber>> graph;
	mutex myMutex,mutexEraseReads;
	atomic<size_t> readMapped,aligneOnPathSucess,unitigsMapped,bigUnitig;



	MappingSupervisor(const vector<string>& Iunitigs, unordered_map<minimizer, unordered_set<rNumber>>& Iindex, size_t Ik, const vector<string>& Ireads, size_t Imulti, size_t IH, size_t Ipart, size_t Ik2, double IminJacc, const unordered_map<string, vector<uNumber>>& Igraph, size_t Ikgraph){
		unitigs=Iunitigs;
		min2Reads=Iindex;
		k=Ik;
		reads=Ireads;
		multi=Imulti;
		H=IH;
		part=Ipart;
		k2=Ik2;
		kgraph=Ikgraph;
		minJacc=IminJacc;
		graph=Igraph;
		readMapped=0;
		aligneOnPathSucess=0;
		unitigsMapped=0;
		offset=100;
		minSizeUnitigs=100;
	}



	void MapPart(size_t L, size_t R);
	void findCandidate(const string& unitig, unordered_set<minimizer>& min, unordered_map<rNumber,size_t>& count, unordered_map<rNumber,unordered_set<minimizer>>& read2min);
	void MapAll();
	int isCandidateCorrect(const string& unitig, rNumber readNumber, unordered_map<rNumber,unordered_set<minimizer>>& read2min, unordered_set<minimizer>& genomicKmers);
	vector<path> listPath(size_t lengthRequired, uNumber ind, unordered_set<uNumber>& usedUnitigs);
	bool alignOnPath(const path& path, const string& read, size_t position,unordered_set<uNumber>& usedUnitigs);
};



#endif /* defined(__PBMOG__MappingSupervisor__) */
