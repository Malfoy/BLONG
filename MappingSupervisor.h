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
	size_t offset,minSizeUnitigs,nbThreads;
	vector<string> unitigs;
	unordered_map<minimizer, vector<rNumber>> min2Reads;
	vector<rPosition> number2position;
	uint32_t readNumber;
	size_t k,multi,H,part,k2,kgraph;
	char depthMax,nuc;
	double minJacc;
	bool mapPartAllowed,errorInKmers,checking;
	double errorRate;
	atomic<double> globalscore;
	graph G;
	mutex mutexReadReads,mutexEraseReads,mutexunitig,mutexReadUsed;
	atomic<size_t> readMapped,aligneOnPathSucess,unitigsPreMapped,bigUnitig,island,regionmapped,leftmap,rightmap,leftmapFail,rightmapFail,candidate,fail,indice,readInUnitig,deepper,failedCompaction,pathNumber,candidateNumber,pathlength;



	MappingSupervisor(const vector<string>& Iunitigs, unordered_map<minimizer, vector<rNumber>>& Iindex, size_t Ik, const string& IreadsFile, size_t Imulti, size_t IH, size_t Ipart, size_t Ik2, double IminJacc,graph& graphe,size_t Ikgraph,vector<rPosition>& Ivect,uint32_t IreadNumber){
		readNumber=IreadNumber;
		number2position=Ivect;
		nuc=4;
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
		offset=200;
		minSizeUnitigs=200;
		depthMax=10;
		globalscore=bigUnitig=regionmapped=deepper=fail=candidate=leftmap=rightmap=leftmapFail=rightmapFail=readInUnitig=failedCompaction=pathNumber=candidateNumber=pathlength=0;
		nbThreads=8;
		errorRate=30;
		indice=0;
		outFile.open("myout.fa",ofstream::trunc);
	}



	void MapPart();
	string getRead(rNumber pos);
	bool mapUnitig(const string& unitig, const rPosition n, int& position, unordered_set<minimizer>& genomicKmers, string& read);
	void findCandidate(const string& unitig, unordered_map<rNumber,uint32_t>& Candidate, unordered_map<rNumber,unordered_set<minimizer>>& read2Min);
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

	string getPathEnd(const vector<uNumber>& numbers);
	string getPathBegin(const vector<uNumber>& numbers);

	void MapFromUnitigsErrors(const string& unitig);

	bool isCandidateCorrectMap(const string& unitig, const string& read, const unordered_multimap<string,string>& genomicKmers,int& position, unordered_set<minimizer>& setMin,int& positionRead);

	bool alignOnPathsSonsErrors(const vector<path>& path, const string& read, size_t position,vector<uNumber>& numbers);
	bool alignOnPathsSons(const vector<path>& Paths, const string& read, size_t position,vector<uNumber>& numbers);
	bool alignOnPathsSonsErrorsAll(const vector<path>& Paths, const string& read, size_t position,vector<uNumber>& numbers);

	bool preMapUnitig(const string& unitig, string& read,const unordered_multimap<string,string>& genomicKmers,int& position,unordered_set<minimizer>& setMin, int& positionUnitig, bool& stranded, int& position1, int& position2);
	bool mapOnGraph(vector<uNumber>& numberBegin, vector<uNumber>& numberEnd,const string& unitig, const string& read, int position,const  vector<path>& ListFathers,const  vector<path>& ListSons, string& beg, string& end);
	string recoverPath(vector<uNumber>& numberBegin, vector<uNumber>& numberEnd,size_t begsize,size_t endsize,size_t readsize,const string& unitig,bool stranded, int position1, int position2);
};



#endif /* defined(__PBMOG__MappingSupervisor__) */
