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
	//foreach unvisited neighbor
	for (size_t i(0); i<indiceNeigboor.size(); ++i){
		if(usedUnitigs.count(indiceNeigboor[i])==0){
			string neigboor(unitigs[indiceNeigboor[i]]);
			if (neigboor.size()>lengthRequired){
				paths.push_back(path{neigboor,indiceNeigboor[i]});
			}else{
				vector<path> paths2(listPath(lengthRequired-neigboor.size()+kgraph, indiceNeigboor[i],usedUnitigs));
				for (size_t j(0);j<paths2.size();++j){

					string str (compaction(paths2[j].str,unitig,kgraph));
					if(!str.empty()){
//						cout<<paths2[j].str<<" "<<unitig<<endl;
//
//						cout<<"begin1 "<<paths2[j].str.substr(0,kgraph)<<" "<<reversecomplement(paths2[j].str.substr(0,kgraph))<<endl;
//
//						cout<<"begin2 "<<unitig.substr(0,kgraph)<<" "<<reversecomplement(unitig.substr(0,kgraph))<<endl;
//
//						cout<<"end1 "<<paths2[j].str.substr(paths2[j].str.size()-kgraph,kgraph)<<" "<<reversecomplement(paths2[j].str.substr(paths2[j].str.size()-kgraph,kgraph))<<endl;
//
//						cout<<"end2 "<<unitig.substr(unitig.size()-kgraph,kgraph)<<" "<<reversecomplement(unitig.substr(unitig.size()-kgraph,kgraph))<<endl;
//						cout<<"wut"<<endl;
//						cin.get();
						paths2[j].str=str;
					}else{
//						cout<<":("<<endl;
					}

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


bool MappingSupervisor::alignOnPath(const path& path, const string& read, size_t position, unordered_set<uNumber>& usedUnitigsInitial){

	if(read.empty()){return false;}
	if(read.size()<=position){return true;}
	int start ((int)position-(int)path.str.size()*1);

	string region(read.substr(max(0,start),path.str.size()*2));
	if(region.size()<offset){return true;}

	unordered_set<minimizer> genomicKmers=allKmerSet(k2,path.str);

	if(jaccard(k2,read,genomicKmers)>minJacc){

		auto list(listPath(offset, path.lastUnitig, usedUnitigsInitial));
//		cout<<list.size()<<endl;
		for(size_t i(0);i<list.size();++i){
			auto usedUnitigs(usedUnitigsInitial);
			if(alignOnPath(list[i], read, position+path.str.size(), usedUnitigs)){
//				cout<<"sucess"<<endl;
				return true;
			}
		}
//		cout<<"fail"<<endl;
	}


	return false;
}




void MappingSupervisor::MapPart(size_t L, size_t R){
	//foreach unitig (sort of)
	for (size_t i(L); i<R; ++i){

		if(unitigsMapped%100==0){
			cout<<unitigsMapped++<<endl;
		}

		string unitig=unitigs[i];
		if(unitig.size()>minSizeUnitigs){
			bigUnitig++;
			unordered_set<uNumber> usedUnitigsShared;
			vector<path> list(listPath(offset, (uNumber)i, usedUnitigsShared));
			bool done(false);
			unordered_set<minimizer> min;
			unordered_map<rNumber,size_t> Candidate;Candidate.set_empty_key(-1);
			unordered_map<rNumber,unordered_set<minimizer>> read2min;read2min.set_empty_key(-1);
			unordered_set<minimizer> genomicKmers;
			findCandidate(unitig,min,Candidate,read2min);
			//foreach reads that could map on the unitig
			for(auto it=Candidate.begin(); it!=Candidate.end(); ++it){
				if(it->second>=multi){
					int position(isCandidateCorrect(unitig,it->first,read2min,genomicKmers));
					if(position!=-1){
						if(!done){
							done=true;
							unitigsMapped++;
							bool mappedOnGraph(false);
							for(size_t ii(0);ii<list.size() and !mappedOnGraph;++ii){
								unordered_set<uNumber> usedUnitigs(usedUnitigsShared);
								path path(list[ii]);
								string read(reads[it->first]);
								if(alignOnPath(path, read.substr(0,position), 0, usedUnitigs)){
									if(alignOnPath(path, read.substr(position), 0, usedUnitigs)){
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





void MappingSupervisor::MapAll(){
	size_t nbThreads(1);
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





