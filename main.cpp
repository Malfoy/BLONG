#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <functional>
#include <unordered_set>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <atomic>
#include "MappingSupervisor.h"
#include "Utils.h"
#include "graph.h"


using namespace std;


int main(){
	size_t H(100),k(15),part(10),kgraph(30),k2(10),threshold(1);
	bool homo(false);
	srand((int)time(NULL));
	size_t nCycle(0);
	double minjacc(10);
	string fileName("/local/malfoyishere/yeast-pacbio.fa");
	cout<<"minjacc : "<<minjacc<<endl;

	auto start=chrono::system_clock::now();
	readContigsforstats("unitigClean.fa", kgraph, false, false, false);
	//~ readContigsforstats("unitigClean.fa", kgraph, false, false, false);
	for(size_t i(0);i<nCycle;++i){
		readContigsforstats("unitigClean.fa", kgraph, true, true, false);
	}
	auto Unitigs(loadUnitigs("unitigClean.fa",homo));
	graph Graph(Unitigs,kgraph);
	//	auto filter(filterUnitigs(Unitigs,k,H,part));
	auto end1=chrono::system_clock::now();auto waitedFor=end1-start;
	cout<<"Unitig loaded "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;

	vector<rPosition> number2position;
	auto index(indexSeqDisk(fileName,H,k,part,number2position));
	auto end2=chrono::system_clock::now();waitedFor=end2-end1;
	cout<<"Reads indexed "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;

	MappingSupervisor supervisor(Unitigs, index, k,fileName, threshold, H, part, k2, minjacc, Graph, kgraph,number2position,number2position.size());
	supervisor.MapAll();
	auto end3=chrono::system_clock::now();waitedFor=end3-end2;
	cout<<"Read aligned in "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;

	return 0;
}
