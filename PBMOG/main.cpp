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

using namespace std;
//using namespace google;

mutex myMutex,myMutex2;
atomic<size_t> atomicount;
uint32_t nindexfasta(0);



void computeMinHash(size_t H, size_t k, size_t part, const vector<string>& sequences, unordered_set<minimizer>* set, size_t L, size_t R){
	for(size_t i(L); i<R; ++i){
		string sequence(sequences[i]);
		if(sequence.size()>=k){
			vector<minimizer> sketch;
			if(sequence.size()<=H){
				//			if(true){
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



void indexSeqAux(const vector<string>& seqs, size_t H, size_t k, size_t part,  const unordered_set<minimizer>& filter, unordered_map<minimizer,vector<rPosition>>* index, uint32_t L, size_t R){
	for (uint32_t i(L); i<R; ++i){
		string seq=seqs[i];
		vector<minimizer> sketch;
		if(seq.size()<=(uint)H){
			sketch=allHash(k,seq);
		}else{
			sketch=minHashpart2(H,k,seq,part,filter);
		}
		removeDuplicate(sketch);
		myMutex.lock();
		for(uint32_t j(0);j<sketch.size();++j){
			(*index)[sketch[j]].push_back(i);
		}
		myMutex.unlock();
	}
}


void indexSeqAux2(const vector<string>& seqs, size_t H, size_t k, size_t part, unordered_map<minimizer,vector<rPosition>>* index, uint32_t L, size_t R){
	for (uint32_t i(L); i<R; ++i){
		string seq=seqs[i];
		vector<minimizer> sketch;
		if(seq.size()<=(uint)H){
			sketch=allHash(k,seq);
		}else{
			sketch=minHashpart(H, k, seq, part);
		}
		removeDuplicate(sketch);
		myMutex.lock();
		for(uint32_t j(0);j<sketch.size();++j){
			(*index)[sketch[j]].push_back(i);
		}
		myMutex.unlock();
	}
}


void indexSeqAux3( vector<string>* seqs, size_t H, size_t k, size_t part, vector<rPosition>* index, uint32_t L, size_t R){
	string seq;
	vector<minimizer> sketch;
	for (uint32_t i(L); i<R; ++i){
		seq=(*seqs)[i];
		if(seq.size()<=(uint)H){
			sketch=allHash(k,seq);
		}else{
			sketch=minHashpart(H, k, seq, part);
		}
		removeDuplicate(sketch);
		myMutex.lock();
		for(uint32_t j(0);j<sketch.size();++j){
			index[sketch[j]].push_back(i);
		}
		myMutex.unlock();
	}
}


unordered_map<minimizer,vector<rPosition>> indexSeq(const vector<string>& seqs, size_t H, size_t k, size_t part, const unordered_set<minimizer>& filter){
	size_t nbThreads(8);
	vector<thread> threads;
	unordered_map<minimizer,vector<rPosition>> index;
	//	index.set_empty_key(-1);
	vector<size_t> limits = bounds(nbThreads, seqs.size());

	for (size_t i(0); i<nbThreads; ++i){
		threads.push_back(thread(indexSeqAux, seqs, H, k, part, filter, &index, limits[i], limits[i+1]));
	}

	for(auto &t : threads){t.join();}
	return index;
}

ifstream global;


void indexFasta(size_t H, size_t k, size_t part, unordered_map<minimizer,vector<rNumber>>* index, int* readNumber, vector<rPosition>* number2position){
	vector<minimizer> sketch;
	string seq,more;
	rPosition position;
	while (!global.eof()){
		bool take(true);
		myMutex2.lock();
		getline(global, seq);
		position=global.tellg();
		getline(global, seq);
		while(global.peek()!='>' and !global.eof()){
			getline(global,more);
			seq+=more;
		}
		myMutex2.unlock();
//		cout<<seq<<endgetl;
//		for (size_t i(0); i<seq.size(); ++i) {
//			if(seq[i]!='A' and seq[i]!='C' and seq[i]!='G' and seq[i]!='T'){
//				take=false;
//				break;
//			}
//		}
		if (take){
			if(seq.size()<=(uint)H){
				sketch=allHash(k,seq);
			}else{
				sketch=minHashpart(H, k, seq, part);
			}
			removeDuplicate(sketch);

			myMutex.lock();
			(*readNumber)++;
			number2position->push_back(position);
			rNumber n((uint32_t)number2position->size()-1);
			for(uint32_t j(0);j<sketch.size();++j){
				(*index)[sketch[j]].push_back(n);
			}
			myMutex.unlock();
//			if(n>10000){return;}
		}
	}
//	cout<<n<<" reads"<<endl;
}


unordered_map<minimizer,vector<rNumber>> indexSeqDisk(const string& seqs, size_t H, size_t k, size_t part,int& readNumber,vector<rPosition>& number2position){
	size_t nbThreads(4);
	vector<thread> threads;
	unordered_map<minimizer,vector<rNumber>> index;
	global.open(seqs);

	for (size_t i(0); i<nbThreads; ++i){
		threads.push_back(thread(indexFasta, H, k, part, &index,&readNumber,&number2position));
	}

	for(auto &t : threads){t.join();}
	return index;
}


vector<rPosition>* indexSeq3( vector<string>* seqs, size_t H, size_t k, size_t part){
	size_t nbThreads(4);
	vector<thread> threads;
	vector<rPosition>* index=new vector<rPosition> [4194304];
	vector<size_t> limits = bounds(nbThreads, seqs->size());

	for (size_t i(0); i<nbThreads; ++i){
		threads.push_back(thread(indexSeqAux3, seqs, H, k, part, index, limits[i], limits[i+1]));
	}

	for(auto &t : threads){t.join();}
	return index;
}






int main(){
	size_t H(100),k(15),part(1),kgraph(30),k2(11),threshold(2);
	bool homo(false);
	srand((int)time(NULL));
	size_t nCycle(0);
	double minjacc(10);
	int readNumber(0);
	string fileName("/Applications/PBMOG/Build/Products/Debug/raw_pacbio_c_elegans.fasta");
	cout<<"minjacc : "<<minjacc<<endl;

	auto start=chrono::system_clock::now();
	readContigsforstats("/Applications/PBMOG/Build/Products/Debug/unitigClean.fa", kgraph, false, false, false);
	for(size_t i(0);i<nCycle;++i){
		readContigsforstats("/Applications/PBMOG/Build/Products/Debug/unitigClean.fa", kgraph, true, true, false);
	}
	auto Unitigs(loadUnitigs("/Applications/PBMOG/Build/Products/Debug/unitigClean.fa",homo));
	graph Graph(Unitigs,kgraph);
	//	auto filter(filterUnitigs(Unitigs,k,H,part));
	auto end1=chrono::system_clock::now();auto waitedFor=end1-start;
	cout<<"Unitig loaded "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;

	//	unordered_map<minimizer, vector<rNumber>> index;
	//	auto vect(indexSeq3(&Reads,H,k,part));
	vector<rPosition> number2position;
	auto index(indexSeqDisk(fileName,H,k,part,readNumber,number2position));
	auto end2=chrono::system_clock::now();waitedFor=end2-end1;
	cout<<"Reads indexed "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;

	MappingSupervisor supervisor(Unitigs, index, k,fileName, threshold, H, part, k2, minjacc, Graph, kgraph,number2position,readNumber);
	supervisor.MapAll();
	auto end3=chrono::system_clock::now();waitedFor=end3-end2;
	cout<<"Read aligned in "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;
	
	return 0;
	
}
