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
	ifstream reads;
	vector<string> unitigs;
	unordered_set<uint> readNumberrecruited;
	unordered_map<minimizer, vector<rNumber>> min2Reads;
	vector<rPosition> number2position;
	uint32_t readNumber;
	unsigned int k,multi,H,part,k2,kgraph,depthMax,nuc,offset,minSizeUnitigs,nbThreads,minJacc;
	bool mapPartAllowed,errorInKmers,checking,first;
	double errorRate;
	atomic<double> globalscore;
	graph G;
	mutex mutexReadReads,mutexEraseReads,mutexunitig,mutexReadUsed;
	atomic<unsigned int> readMapped,aligneOnPathSucess,unitigsPreMapped,bigUnitig,island,regionmapped,leftmap,rightmap,leftmapFail,rightmapFail,candidate,fail,indice,readInUnitig,deepper,failedCompaction,pathNumber,candidateNumber,pathlength;



	MappingSupervisor(const vector<string>& Iunitigs, unordered_map<minimizer, vector<rNumber>>& Iindex, uint  Ik, const string& IreadsFile, uint  Imulti, uint  IH, uint  Ipart, uint  Ik2, double IminJacc,graph& graphe,uint  Ikgraph,vector<rPosition>& Ivect,uint32_t IreadNumber,uint  Ithread,const string& Ioutput, uint  Ioffset,uint  Iminsizeunitig,uint  Idepth){
		readNumber=IreadNumber;
		number2position=Ivect;
		nuc=5;
		readNumberrecruited={};
		mapPartAllowed=false;
		errorInKmers=true;
		checking=false;
		unitigs=Iunitigs;
		min2Reads=Iindex;
		k=Ik;
		reads.open(IreadsFile);
//		reads=Ireads;
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
		offset=Ioffset;
		minSizeUnitigs=Iminsizeunitig;
		depthMax=Idepth;
		globalscore=bigUnitig=regionmapped=deepper=fail=candidate=leftmap=rightmap=leftmapFail=rightmapFail=readInUnitig=failedCompaction=pathNumber=candidateNumber=pathlength=0;
		nbThreads=Ithread;
		errorRate=10;
		indice=0;
		first=true;
		outFile.open(Ioutput,ofstream::trunc);
	}


	void MapPart();
	string getRead(rNumber pos);
	bool mapUnitig(const string& unitig, const rPosition n, int& position, unordered_set<minimizer>& genomicKmers, string& read);
	void findCandidate(const string& unitig, unordered_map<rNumber,uint32_t>& Candidate);
	void MapAll();
	bool isCandidateCorrect(const string& unitig, const string& read, unordered_set<minimizer>& genomicKmers,int& position, unordered_set<minimizer>& setmin);
	bool isCandidateCorrect2(const string& unitig, const string& read, unordered_set<minimizer>& genomicKmers,int& position, unordered_set<minimizer>& setmin);
	vector<path> listPath(uint  lengthRequired, uNumber ind, set<uNumber> usedUnitigs);

	bool alignOnPath(const path& path, const string& read, uint  position,set<uNumber>& usedUnitigs);

	vector<path> listPathSons(uint  lengthRequired, const string& substr,uint depth);
	vector<path> listPathFathers(uint  lengthRequired, const string& substr, uint depth);

	bool alignOnPathSons(const path& path, const string& read, uint  position,vector<uNumber>& n);
	bool alignOnPathSonsErrors(const path& path, const string& read, uint  position,vector<uNumber>& n);
	bool alignOnPathFathers(const path& path, const string& read, uint  position,vector<uNumber>& n);

	bool alignOnPathSons2(const path& path, const string& read, uint  position,vector<uNumber>&);
	bool alignOnPathFathers2(const path& path, const string& read, uint  position,vector<uNumber>&);

	string getPathEnd(const vector<uNumber>& numbers);
	string getPathBegin(const vector<uNumber>& numbers);

	void MapFromUnitigsErrors(const string& unitig);

	bool isCandidateCorrectMap(const string& unitig, const string& read, vector<uint32_t>* genomicKmers,int& position, unordered_set<minimizer>& setMin, int& positionRead);
	bool alignOnPathsSonsErrors(const vector<path>& path, const string& read, uint  position,vector<uNumber>& numbers);
	bool alignOnPathsSons(const vector<path>& Paths, const string& read, uint  position,vector<uNumber>& numbers);
	bool alignOnPathsSonsErrorsAll(const vector<path>& Paths, const string& read, uint  position,vector<uNumber>& numbers);

	bool preMapUnitig(const string& unitig, string& read,vector<uint32_t>* genomicKmers,int& position,unordered_set<minimizer>& setMin, int& positionUnitig, bool& stranded, int& position1);
	bool mapOnGraph(vector<uNumber>& numberBegin, vector<uNumber>& numberEnd,const string& unitig, const string& read, int position,const  vector<path>& ListFathers,const  vector<path>& ListSons, string& beg, string& end);
	string recoverPath(vector<uNumber>& numberBegin, vector<uNumber>& numberEnd,const string& unitig);

	unordered_set<minimizer> allKmerset(uint  k,const string& seq);

};



#endif /* defined(__PBMOG__MappingSupervisor__) */
