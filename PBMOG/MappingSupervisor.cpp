//
//  MappingSupervisor.cpp
//  PBMOG
//
//  Created by malfoy on 09/01/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#include "MappingSupervisor.h"
#include "Utils.h"

void MappingSupervisor::MapPart(size_t L, size_t R){
	size_t minSizeUnitigs(100);

	for (size_t i(L); i<R; ++i){
		string unitig=unitigs[i];
		if(unitig.size()>minSizeUnitigs){
			unordered_set<minimizer> min;
			unordered_map<rNumber,size_t> count;
			unordered_map<rNumber,unordered_set<minimizer>> read2min;
			unordered_set<minimizer> A;
			if(unitig.size()<=H){
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
				vector<minimizer> sketch;
				sketch=minHashpart(H,k,unitig,part);

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
			for(auto it=count.begin();it!=count.end();++it){
				uint32_t readNumber=it->first;
				bool goodreadb=false;
				if(it->second>=multi){
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
							if(A.size()==0){
								vector<minimizer> sketchUnitig=allHash(k2,unitig);
								for(size_t m(0);m<sketchUnitig.size();++m){
									A.insert(sketchUnitig[m]);
								}
							}
							if(jaccard3(k2,region,A)>minJacc){
								atomicount++;
							}
						}
					}
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