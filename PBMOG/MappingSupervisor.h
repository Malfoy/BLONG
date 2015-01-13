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

class MappingSupervisor{
public:
	vector<string> unitigs,reads;
	unordered_map<minimizer, unordered_set<rNumber>> index;
	size_t k,multi,H,part,k2;
	double minJacc;
	unordered_map<string, vector<uint32_t>> graph;
	mutex myMutex;
	atomic<size_t> atomicount;


	MappingSupervisor(const vector<string>& Iunitigs, unordered_map<minimizer, unordered_set<rNumber>>& Iindex, size_t Ik, const vector<string>& Ireads, size_t Imulti, size_t IH, size_t Ipart, size_t Ik2, double IminJacc, const unordered_map<string, vector<uint32_t>>& Igraph){
		unitigs=Iunitigs;
		index=Iindex;
		k=Ik;
		reads=Ireads;
		multi=Imulti;
		H=IH;
		part=Ipart;
		k2=Ik2;
		minJacc=IminJacc;
		graph=Igraph;
	}

	void MapPart(size_t L, size_t R);
};



#endif /* defined(__PBMOG__MappingSupervisor__) */
