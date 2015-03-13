#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <functional>
#include <unordered_set>
//#include <sparsehash/sparse_hash_map>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <atomic>
#include "MappingSupervisor.h"
#include "Utils.h"
#include "binSeq.h"
#include "graph.h"

//~ #define unordered_map sparse_hash_map
#define minimizer uint32_t
#define rNumber uint32_t

using namespace std;
//using namespace google;

mutex myMutex;
atomic<size_t> atomicount;



void computeMinHash(size_t H, size_t k, size_t part, const vector<string>& sequences, unordered_set<minimizer>* set, size_t L, size_t R){
	for(size_t i(L); i<R; ++i){
		string sequence(sequences[i]);
		if(sequence.size()>=k){
			vector<minimizer> sketch;
			if(sequence.size()<=H){
				sketch=allHash(k, sequence);
			}else{
				sketch=minHashpart(H,k,sequence,part);
			}
			myMutex.lock();
			for(size_t j(0);j<sketch.size();++j){
				set->insert(sketch[j]);
			}
			myMutex.unlock();
		}
	}
}




unordered_set<minimizer> filterUnitigs(const vector<string>& unitigs, size_t k, size_t H, size_t part){
	size_t nbThreads(8);
	unordered_set<minimizer> res;
	string line;
	vector<thread> threads;

	vector<size_t> limits = bounds(nbThreads, unitigs.size());

	for (size_t i(0); i<nbThreads; ++i){
		threads.push_back(thread(computeMinHash,H,k,part,unitigs,&res,limits[i],limits[i+1]));
	}

	for(auto &t : threads){t.join();}
	return res;
}



void indexSeqAux(const vector<string>& seqs, size_t H, size_t k, size_t part,  const unordered_set<minimizer>& filter, unordered_map<minimizer,vector<rNumber>>* index, uint32_t L, size_t R){

	for (uint32_t i(L); i<R; ++i){
		string seq=seqs[i];
		vector<minimizer> sketch;
		if(seq.size()<=(uint)H){
			sketch=allHash(k,seq);
		}else{
			sketch=minHashpart2(H,k,seq,part,filter);
		}
		myMutex.lock();
		for(uint32_t j(0);j<sketch.size();++j){
			(*index)[sketch[j]].push_back(i);
		}
		myMutex.unlock();
	}
}



unordered_map<minimizer,vector<rNumber>> indexSeq(const vector<string>& seqs, size_t H, size_t k, size_t part, const unordered_set<minimizer>& filter){
	size_t nbThreads(8);
	vector<thread> threads;
	unordered_map<minimizer,vector<rNumber>> index;
	index.set_empty_key(-1);
	vector<size_t> limits = bounds(nbThreads, seqs.size());

	for (size_t i(0); i<nbThreads; ++i) {
		threads.push_back(thread(indexSeqAux, seqs, H, k, part, filter, &index, limits[i], limits[i+1]));
	}

	for(auto &t : threads){t.join();}
	return index;
}



void testSimilarity(const string& refFaFile, const string& pbFileFq){
	size_t k(11);
	bool homo(false);
	ifstream refFile(refFaFile),pbFile(pbFileFq);
	string read,ref,line;
	getline(refFile,line);
	getline(refFile,ref);

	vector<string> pbReads;

	while (!pbFile.eof()) {
		getline(pbFile,line);
		getline(pbFile,read);
		getline(pbFile,line);
		getline(pbFile,line);
		pbReads.push_back(read);
	}
	unordered_set<minimizer> genomicKmers;
	if(homo){
		genomicKmers=allKmerSet(k,homocompression(ref));
	}else{
		genomicKmers=allKmerSet(k,ref);
	}
	//	unordered_set<minimizer> genomicKmers=allKmerSet(k2,path.str);


	double acc(0);
	for (size_t i(0); i<pbReads.size(); ++i) {
		read=pbReads[i];
		if(homo){
			double jacc(jaccard(k,homocompression(read),genomicKmers));
			acc+=jacc;
			cout<<jacc<<endl;
		}else{
			double jacc(jaccard(k,read,genomicKmers));
			cout<<jacc<<endl;
			acc+=jacc;
		}
	}
	cout<<"read number : "<<pbReads.size()<<endl;
	cout<<"mean : "<<acc/pbReads.size()<<endl;
}



int main(){
	//	testBinSeq();
	//	exit(0);
	//	testSimilarity("/Applications/PBMOG/Build/Products/Debug/random.fa","/Applications/PBMOG/Build/Products/Debug/sd_0001.fastq");
	//	exit(0);

//	size_t H(100),k(15),part(1),kgraph(30),k2(7),minsize(100),threshold(3);
	size_t H(100),k(5),part(1),kgraph(5),k2(5),minsize(1),threshold(1);
	bool homo(false);
	srand((int)time(NULL));
	size_t nCycle(0);
	double errorRate(0);
	double minjacc(100*(pow(1-errorRate,k2)));
	//	double minjacc(20);
	cout<<"minjacc : "<<minjacc<<endl;

	auto start=chrono::system_clock::now();
//	auto Reads(loadFASTQ("/Applications/PBMOG/Build/Products/Debug/positive_0001.fastq",homo,minsize));
//	auto Reads(loadFASTQ("/Applications/PBMOG/Build/Products/Debug/1Xnormal10K09_0001.fastq",homo,minsize));
//	readContigsforstats("/Applications/PBMOG/Build/Products/Debug/unitig31.fa", kgraph, false, true, false);
	auto Reads(loadFASTQ("/Applications/PBMOG/Build/Products/Debug/read.fa",homo,minsize));
	readContigsforstats("/Applications/PBMOG/Build/Products/Debug/unitigs.fa", kgraph, false, false, false);
	for(size_t i(0);i<nCycle;++i){
		readContigsforstats("/Applications/PBMOG/Build/Products/Debug/unitigClean.fa", kgraph, true, true, false);
	}

	auto Unitigs(loadUnitigs("/Applications/PBMOG/Build/Products/Debug/unitigClean.fa",homo));
//	vector<string> Unitigs(kmerCounting("/Applications/PBMOG/Build/Products/Debug/ecoliref.fasta", kgraph+1));
	//	auto Unitigs(loadUnitigs("/Applications/PBMOG/Build/Products/Debug/randomref1.fa",homo));
	graph Graph(Unitigs,kgraph);
	auto Filter(filterUnitigs(Unitigs,k,H,part));
	auto end1=chrono::system_clock::now();auto waitedFor=end1-start;
	cout<<"Reads loaded "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;

	auto index(indexSeq(Reads,H,k,1*part,Filter));
	auto end2=chrono::system_clock::now();waitedFor=end2-end1;
	cout<<"Reads indexed "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;

	MappingSupervisor supervisor(Unitigs, index, k, Reads, threshold, H, part, k2, minjacc, Graph, kgraph);

	supervisor.MapAll();
	auto end3=chrono::system_clock::now();waitedFor=end3-end2;
	cout<<"Reads Mapped "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;

	return 0;
}
