//
//  MappingSupervisor.cpp
//  PBMOG
//
//  Created by malfoy on 09/01/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#include "MappingSupervisor.h"
#include "Utils.h"





vector<path> MappingSupervisor::listPath(size_t lengthRequired, uNumber ind, unordered_set<uNumber>& usedUnitigs){
	string unitig(unitigs[ind]);
	vector<path> paths;
	vector<uNumber> indiceNeigboor(graph[getRepresent(unitig.substr(0,kgraph))]);
	vector<uNumber> indiceNeigboor2(graph[getRepresent(unitig.substr(unitig.size()-kgraph,kgraph))]);
	indiceNeigboor.insert(indiceNeigboor.end(), indiceNeigboor2.begin(), indiceNeigboor2.end());
	usedUnitigs.insert(ind);
	for (size_t i(0); i<indiceNeigboor.size(); ++i){
		if(usedUnitigs.count(indiceNeigboor[i])==0){
			string neigboor(unitigs[indiceNeigboor[i]]);
			if (neigboor.size()>lengthRequired){
				paths.push_back(path{neigboor,indiceNeigboor[i]});
			}else{
				vector<path> paths2(listPath(lengthRequired-neigboor.size()+kgraph, indiceNeigboor[i],usedUnitigs));
				for (size_t j(0);j<paths2.size();++j){
					paths2[j].str=compaction(paths2[j].str,unitig,kgraph);
				}
				paths.insert(paths.end(), paths2.begin(), paths2.end());
			}
		}
	}
	return paths;
}




void MappingSupervisor::findCandidate(const string& unitig, unordered_set<minimizer>& min, unordered_map<rNumber,size_t>& Candidate, unordered_map<rNumber,unordered_set<minimizer>>& read2Min){
	if(unitigs.size()<100){
		for(size_t j(0);j+k<unitig.size();++j){
			minimizer seq=seq2int(unitig.substr(j,k));
			if(min.count(seq)==0){
				if(min2Reads.count(seq)!=0){
					min.insert(seq);
					unordered_set<rNumber> myset=min2Reads[seq];
					for(auto it=myset.begin();it!=myset.end();++it){
						if(!reads[*it].empty()){
							Candidate[*it]++;
							read2Min[*it].insert(seq);
						}
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
				if(min2Reads.count(seq)!=0){
					unordered_set<rNumber> myset=min2Reads[seq];
					for(auto it=myset.begin();it!=myset.end();++it){
						if(!reads[*it].empty()){
							++Candidate[*it];
							read2Min[*it].insert(seq);
						}
					}
				}
			}
		}
	}
}




int MappingSupervisor::isCandidateCorrect(const string& unitig, rNumber readNumber, unordered_map<rNumber,unordered_set<minimizer>>& read2min, unordered_set<minimizer>& genomicKmers){
	bool goodreadb=false;
	unordered_set<minimizer> setMin(read2min[readNumber]);
	string read=reads[readNumber];
	if(read.empty()){return -1;}
	for(size_t j(0);j+k<read.size() and !goodreadb;++j){
		minimizer seq(seq2int(read.substr(j,k)));
		if(setMin.count(seq)!=0){
			string region;

			if(j>1*unitig.size()){
				region=read.substr(j-1*unitig.size(),2*unitig.size());
			}else{
				region=read.substr(0,2*unitig.size());
			}
			if(genomicKmers.size()==0){
				genomicKmers=allKmerSet(k2,unitig);
			}
			if(jaccard(k2,region,genomicKmers)>minJacc){
				readMapped++;
				return ((int)j);
			}
		}
	}
	return -1;
}


bool MappingSupervisor::alignOnPath(const path& path, const string& read, size_t position, unordered_set<uNumber>& usedUnitigs){

	if(read.empty()){return false;}
	if(read.size()<position){return true;}

	string region(read.substr(position,path.str.size()));
	if(region.size()<100){return true;}

	unordered_set<minimizer> genomicKmers=allKmerSet(k2,path.str);
//	if(jaccard(k2,read,genomicKmers)>minJacc){
//		aligneOnPathSucess++;
//		return true;
//	}

	if(jaccard(k2,read,genomicKmers)>minJacc){
		auto list(listPath(100, path.lastUnitig, usedUnitigs));
		for(size_t i(0);i<list.size();++i){
			if(alignOnPath(list[i], read, position+100, usedUnitigs)){
				return true;
			}
		}
	}


	return false;
}




void MappingSupervisor::MapPart(size_t L, size_t R){
	size_t minSizeUnitigs(100);
	for (size_t i(L); i<R; ++i){
		if(unitigsMapped%100==0){
			cout<<unitigsMapped++<<endl;
		}
		string unitig=unitigs[i];
		if(unitig.size()>minSizeUnitigs){
			bigUnitig++;
			unordered_set<uNumber> usedUnitigs;
			vector<path> list(listPath(1*100, (uNumber)i, usedUnitigs));
			if(unitig.size()>minSizeUnitigs){
				bool done(false);
				unordered_set<minimizer> min;
				unordered_map<rNumber,size_t> Candidate;
				Candidate.set_empty_key(-1);
				unordered_map<rNumber,unordered_set<minimizer>> read2min;
				read2min.set_empty_key(-1);
				unordered_set<minimizer> genomicKmers;
				findCandidate(unitig,min,Candidate,read2min);
				for(auto it=Candidate.begin(); it!=Candidate.end(); ++it){
					if(it->second>=multi){
						int position(isCandidateCorrect(unitig,it->first,read2min,genomicKmers));
						if(position!=-1){
							if(!done){
								done=true;
								unitigsMapped++;
								bool mappedOnGraph(false);
								for(size_t ii(0);ii<list.size() and !mappedOnGraph;++ii){
									path path(list[ii]);
									if(alignOnPath(path, reads[it->first].substr(0,position),0,usedUnitigs)){
										if(alignOnPath(path, reads[it->first].substr(position), 0,usedUnitigs)){
											aligneOnPathSucess++;
											mappedOnGraph=true;
											mutexEraseReads.lock();
											reads[it->first].clear();
											mutexEraseReads.unlock();
										}
									}
								}
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
	vector<size_t> limits = bounds(nbThreads, unitigs.size());

	for(size_t i(0); i<nbThreads; ++i) {
		threads.push_back(thread(&MappingSupervisor::MapPart, this, limits[i], limits[i+1]));
	}

	for(auto &t : threads){t.join();}
	cout<<"Unitigs mapped "<<unitigsMapped<<" Percent unitigs mapped : "<<(100*unitigsMapped)/(bigUnitig)<<endl;
	cout<<"Read mapped : "<<readMapped<<" Percent read mapped : "<<(100*readMapped)/(reads.size()+1)<<endl;
	cout<<"Read mapped on graph: "<<aligneOnPathSucess<<" Percent read mapped : "<<(100*aligneOnPathSucess)/(reads.size()+1)<<endl;
}





