//
//  MappingSupervisor.cpp
//  PBMOG
//
//  Created by malfoy on 09/01/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#include "MappingSupervisor.h"
#include "Utils.h"

void MappingSupervisor::findCandidate(const string& unitig, unordered_set<minimizer>& min, unordered_map<rNumber,size_t>& count, unordered_map<rNumber,unordered_set<minimizer>>& read2min){
	if(unitigs.size()<100){
		for(size_t j(0);j+k<unitig.size();++j){
			uint32_t seq=seq2int(unitig.substr(j,k));
			if(min.count(seq)==0){
				if(index.count(seq)!=0){
					min.insert(seq);
					unordered_set<rNumber> myset=index[seq];
					for(auto it=myset.begin();it!=myset.end();++it){
						count[*it]++;
						read2min[*it].insert(seq);
					}
				}
			}
		}
	}else{
		vector<minimizer> sketch(minHashpart(H,k,unitig,part));
		//For each minimizer of the unitig
		for(unsigned int j(0);j<sketch.size();++j){
			minimizer seq=sketch[j];
			if(min.count(seq)==0){
				min.insert(seq);
				if(index.count(seq)!=0){
					unordered_set<rNumber> myset=index[seq];
					for(auto it=myset.begin();it!=myset.end();++it){
						++count[*it];
						read2min[*it].insert(seq);
					}
				}
			}
		}
	}
}

void MappingSupervisor::isCandidateCorrect(const string& unitig, uint32_t readNumber,unordered_map<rNumber,unordered_set<minimizer>> read2min,unordered_set<minimizer> genomicKmers){
	bool goodreadb=false;
	unordered_set<minimizer> setMin(read2min[readNumber]);
	string read=reads[readNumber];
	for(uint j(0);j+k<read.size() and !goodreadb;++j){
		minimizer seq(seq2int(read.substr(j,k)));
		if(setMin.count(seq)!=0){
			string region;

			if(j>1*unitig.size()){
				region=read.substr(j-1*unitig.size(),2*unitig.size());
			}else{
				region=read.substr(0,2*unitig.size());
			}
			if(genomicKmers.size()==0){
				vector<minimizer> sketchUnitig=allHash(k2,unitig);
				for(size_t m(0);m<sketchUnitig.size();++m){
					genomicKmers.insert(sketchUnitig[m]);
				}
			}
			if(jaccard(k2,region,genomicKmers)>minJacc){
				atomicount++;
				goodreadb=true;
			}
		}
	}
}


void MappingSupervisor::MapPart(size_t L, size_t R){
	size_t minSizeUnitigs(100);

	for (size_t i(L); i<R; ++i){
		string unitig=unitigs[i];
		if(unitig.size()>minSizeUnitigs){
			unordered_set<minimizer> min;
			unordered_map<rNumber,size_t> Candidate;
			unordered_map<rNumber,unordered_set<minimizer>> read2min;
			unordered_set<minimizer> genomicKmers;
			findCandidate(unitig,min,Candidate,read2min);

			for(auto it=Candidate.begin();it!=Candidate.end();++it){
				if(it->second>=multi){
					isCandidateCorrect(unitig,it->first,read2min,genomicKmers);
				}
			}
		}
	}
}




void MappingSupervisor::MapAll(){
	size_t nbThreads(8);
	vector<thread> threads;
	vector<size_t> limits = bounds(nbThreads, reads.size());

	for (size_t i(0); i<nbThreads; ++i) {
		threads.push_back(thread(&MappingSupervisor::MapPart, this, limits[i], limits[i+1]));
	}

	for(auto &t : threads){t.join();}

	cout<<"Read mapped : "<<atomicount<<" Percent read mapped : "<<atomicount/(reads.size()+1)<<endl;
}