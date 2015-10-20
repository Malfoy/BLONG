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
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


using namespace std;


int main(int argc, char ** argv){
	size_t H(100),k(15),part(2),kgraph(30),k2(10),threshold(2),nCycle(0),minjacc(10);
	string readName(""),unitigName(""),outFile("out.fa");
	int c;
	while ((c = getopt (argc, argv, "r:u:c:H:K:p:t:k:m:o")) != -1){
		switch(c){
			case 'r':
				readName=optarg;
				break;
			case 'u':
				unitigName=optarg;
				break;
			case 'o':
				outFile=optarg;
				break;
			case 'c':
				nCycle=stoi(optarg);
				break;
			case 'H':
				H=stoi(optarg);
				break;
			case 'K':
				k=stoi(optarg);
				break;
			case 'p':
				part=stoi(optarg);
				break;
			case 't' :
				threshold=stoi(optarg);
				break;
			case 'k' :
				k2=stoi(optarg);
				break;
			case 'm' :
				minjacc=stoi(optarg);
				break;
		}
	}

	if(!unitigName.empty() and !readName.empty()){
		srand((int)time(NULL));
		auto start=chrono::system_clock::now();
		readContigsforstats(unitigName, kgraph, false, false, true);
		for(size_t i(0);i<nCycle;++i){
			readContigsforstats("unitigClean.fa", kgraph, true, true, false);
		}
		auto Unitigs(loadUnitigs("unitigClean.fa",false));
		graph Graph(Unitigs,kgraph);
		//	auto filter(filterUnitigs(Unitigs,k,H,part));
		auto end1=chrono::system_clock::now();auto waitedFor=end1-start;
		cout<<"Unitig loaded "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;

		vector<rPosition> number2position;
		auto index(indexSeqDisk(readName,H,k,part,number2position));
		auto end2=chrono::system_clock::now();waitedFor=end2-end1;
		cout<<"Reads indexed "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;

		MappingSupervisor supervisor(Unitigs, index, k,readName, threshold, H, part, k2, minjacc, Graph, kgraph,number2position,number2position.size());
		supervisor.MapAll();
		auto end3=chrono::system_clock::now();waitedFor=end3-end2;
		cout<<"Read aligned in "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;
	}else{
		cout<<"provide read and unitig file ..."<<endl;
	}

	return 0;
}
