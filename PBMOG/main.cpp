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

//~ #define unordered_map sparse_hash_map
#define minimizer uint32_t
#define rNumber uint32_t

using namespace std;
//using namespace google;

mutex myMutex;
atomic<size_t> atomicount;



void computeMinHash(size_t H, size_t k, size_t part, const vector<string>& reads, unordered_set<minimizer>* set, size_t L, size_t R){
	for(size_t i(L); i<R; ++i){
		string read(reads[i]);
		if(read.size()>100){
			vector<minimizer> sketch=minHashpart(H,k,read,part);
			myMutex.lock();
			for(size_t j(0);j<H;++j){
				set->insert(sketch[j]);
			}
			myMutex.unlock();
		}
	}
}




unordered_set<minimizer> filterUnitigs(const vector<string>& V, size_t k, size_t H, size_t part){
	size_t nbThreads(8);
	unordered_set<minimizer> res;
	string line;
	vector<thread> threads;

	vector<size_t> limits = bounds(nbThreads, V.size());

	for (size_t i(0); i<nbThreads; ++i){
		threads.push_back(thread(computeMinHash,H,k,part,V,&res,limits[i],limits[i+1]));
	}

	for(auto &t : threads){t.join();}
	return res;
}



void indexSeqAux(const vector<string>& seqs, size_t H, size_t k, size_t part,  const unordered_set<minimizer>& filter, unordered_map<minimizer,unordered_set<rNumber>>* index, uint32_t L, size_t R){

	for (uint32_t i(L); i<R; ++i){
		string seq=seqs[i];
		vector<minimizer> sketch;
		if(seq.size()<=(uint)H){
			sketch=allHash(k,seq);
		}else{
			sketch=minHashpart2(H,k,seq,part,filter);
			//~ sketch=minHash(H,k,unitig);
		}
		myMutex.lock();
		for(uint32_t j(0);j<sketch.size();++j){
			(*index)[sketch[j]].insert(i);
		}
		myMutex.unlock();
	}
}



unordered_map<minimizer,unordered_set<rNumber>> indexSeq(const vector<string>& seqs, size_t H, size_t k, size_t part, const unordered_set<minimizer>& filter){
	size_t nbThreads(8);
	vector<thread> threads;
	unordered_map<minimizer,unordered_set<rNumber>> index;
	index.set_empty_key(-1);
	vector<size_t> limits = bounds(nbThreads, seqs.size());

	for (size_t i(0); i<nbThreads; ++i) {
		threads.push_back(thread(indexSeqAux, seqs, H, k, part, filter, &index, limits[i], limits[i+1]));
	}

	for(auto &t : threads){t.join();}
	return index;
}



int main(int argc, char ** argv) {
	size_t H1(100),k(15),part(1),kgraph(19);
	size_t k2(11);
	bool homo(false);
	srand((int)time(NULL));
	size_t nCycle(10);
	double errorRate(0.15);

	auto start=chrono::system_clock::now();
	auto R(loadFASTQ("/Applications/PBMOG/Build/Products/Debug/PC10x_0001.fastq",homo));
	readContigsforstats("/Applications/PBMOG/Build/Products/Debug/unitigEcoliK20NKS2.fa", kgraph, false, true, false);
	for(size_t i(0);i<nCycle;++i){
		readContigsforstats("/Applications/PBMOG/Build/Products/Debug/unitigClean.fa", kgraph, true, true, false);
	}

	auto U(loadUnitigs("/Applications/PBMOG/Build/Products/Debug/unitigClean.fa",homo));
	auto G(getGraph(U,kgraph));
	auto F(filterUnitigs(U,k,H1,part));
	auto end1=chrono::system_clock::now();auto waitedFor=end1-start;
	cout<<"Reads loaded "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;

	auto index(indexSeq(R,H1,k,1*part,F));
	auto end2=chrono::system_clock::now();waitedFor=end2-end1;
	cout<<"Reads indexed "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;

	MappingSupervisor supervisor(U, index, k, R, 2, H1, part, k2, 100*(pow(1-2*errorRate,k2)), G, kgraph);
	supervisor.MapAll();
	auto end3=chrono::system_clock::now();waitedFor=end3-end2;
	cout<<"Reads Mapped "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;

	return 0;
}
