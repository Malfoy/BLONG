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
	uint32_t lastUnitig;
};

class MappingSupervisor{
public:
	vector<string> unitigs,reads;
	unordered_map<minimizer, unordered_set<rNumber>> min2Reads;
	size_t k,multi,H,part,k2,kgraph;
	double minJacc;
	unordered_map<string, vector<uint32_t>> graph;
	mutex myMutex,mutexEraseReads;
	atomic<size_t> readMapped,aligneOnPathSucess,unitigsMapped,bigUnitig;



	MappingSupervisor(const vector<string>& Iunitigs, unordered_map<minimizer, unordered_set<rNumber>>& Iindex, size_t Ik, const vector<string>& Ireads, size_t Imulti, size_t IH, size_t Ipart, size_t Ik2, double IminJacc, const unordered_map<string, vector<uint32_t>>& Igraph, size_t Ikgraph){
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
	}



	void MapPart(size_t L, size_t R);
	void findCandidate(const string& unitig, unordered_set<minimizer>& min, unordered_map<rNumber,size_t>& count, unordered_map<rNumber,unordered_set<minimizer>>& read2min);
	void MapAll();
	int isCandidateCorrect(const string& unitig, uint32_t readNumber, unordered_map<rNumber,unordered_set<minimizer>>& read2min, unordered_set<minimizer>& genomicKmers, uint32_t ind);
	vector<path> listPath(size_t lengthRequired, uint32_t ind, unordered_set<uint32_t>& usedUnitigs);
	bool alignOnPath(const path& path, const string& read, size_t position,unordered_set<uint32_t>& usedUnitigs);
};



#endif /* defined(__PBMOG__MappingSupervisor__) */
