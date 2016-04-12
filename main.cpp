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
	uint  H(100),k(15),part(1),kgraph(62),k2(10),threshold(2),nCycle(0),minjacc(10),thread(1),minSizeUnitig(100),offset(100),depth(10);
	string readName(""),unitigName(""),outFile("out.fa");
	int c;
	while ((c = getopt (argc, argv, "r:u:c:H:K:p:t:k:m:o:t:s:f:d:")) != -1){
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
			case 'n' :
				threshold=stoi(optarg);
				break;
			case 'k' :
				k2=stoi(optarg);
				break;
			case 't' :
				thread=stoi(optarg);
				break;
			case 'm' :
				minjacc=stoi(optarg);
				break;
			case 's' :
				minSizeUnitig=stoi(optarg);
				break;
			case 'f' :
				offset=stoi(optarg);
				break;
			case 'd' :
				depth=stoi(optarg);
				break;
		}
	}

	if(!unitigName.empty() and !readName.empty()){
		cout<<"rReadFile:"<<readName<<"UnitigFile:"<<unitigName<<"cycleNumber:"<<nCycle<<"H:"<<H<<"K:"<<k<<"unitigSize:"<<minSizeUnitig<<"minjacc:"<<minjacc<<"offset:"<<offset<<"smallk:"<<k2<<"depthmax:"<<depth<<endl;
		srand((int)time(NULL));
		auto start=chrono::system_clock::now();
		readContigsforstats(unitigName, kgraph, false, false, false);
		for(uint  i(0);i<nCycle;++i){
			readContigsforstats("data/unitigClean.fa", kgraph, true, true, false);
		}
		cout<<"contigsread"<<endl;
		auto Unitigs(loadUnitigs("data/unitigClean.fa",false));
		graph Graph(Unitigs,kgraph);
		auto filter(filterUnitigs(Unitigs,k,H,part));

		auto end1=chrono::system_clock::now();auto waitedFor=end1-start;
		cout<<"Unitig loaded "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;

		vector<rPosition> number2position;
		auto index(indexSeqDisk(readName,H,k,part,number2position,filter));
		auto end2=chrono::system_clock::now();waitedFor=end2-end1;
		cout<<"Reads indexed "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;

		MappingSupervisor supervisor(Unitigs, index, k,readName, threshold, H, part, k2, minjacc, Graph, kgraph,number2position,number2position.size(),thread,outFile,offset,minSizeUnitig,depth);
		supervisor.MapAll();
		auto end3=chrono::system_clock::now();waitedFor=end3-end2;
		cout<<"Read aligned in "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;
	}else{
		cout<<"Provide at least read and unitig file ..."<<endl;
		cout<<"Option list: "<<endl<<endl;

		cout<<"I/O:"<<endl;
		cout<<"-r Read file"<<endl;
		cout<<"-u Unitig file"<<endl;
		cout<<"-o Out file (out.fa)"<<endl;
		cout<<"-c Number of unitig cleaning operation (0)"<<endl<<endl;

		cout<<"Anchoring:"<<endl;
		cout<<"-H Number minimizer used for anchoring, higher is more sensitive but more expensive (100)"<<endl;
		cout<<"-K Size of  minimizer used for anchoring (15)"<<endl;
		cout<<"-n Minimum number of minimizer shared to be candidate (2)"<<endl;
		cout<<"-s Minimum size for a  unitig to be used as an anchor (100)"<<endl<<endl;;


		cout<<"Matching:"<<endl;
		cout<<"-m Percentage of accepted kmer required for matching (20)"<<endl;
		cout<<"-k Size of kmer used for matching (10)"<<endl;
		cout<<"-f Size of the windows used for the choice of the path (100)"<<endl;
		cout<<"-d Maximal depth in the DFS search (10)"<<endl<<endl;

		cout<<"Performances:"<<endl;
		cout<<"-t Number of thread used (1)"<<endl;
	}

	return 0;
}
