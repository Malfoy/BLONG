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
	vector<uNumber> indiceNeigboor(G.getRight(unitig.substr(0,kgraph)));
	vector<uNumber> indiceNeigboor2(G.getLeft(unitig.substr(unitig.size()-kgraph,kgraph)));
//	cout<<indiceNeigboor.size()<<" "<<indiceNeigboor2.size()<<endl;
//	cout<<getRepresent(unitig.substr(0,kgraph))<<endl;
//	for(size_t i(0);i<indiceNeigboor.size();++i){
//		cout<<unitigs[indiceNeigboor[i]]<<endl;
//
//	}
//	cin.get();
	indiceNeigboor.insert(indiceNeigboor.end(), indiceNeigboor2.begin(), indiceNeigboor2.end());
	usedUnitigs.insert(ind);
	//foreach unvisited neighbor
	for (size_t i(0); i<indiceNeigboor.size(); ++i){
		if(usedUnitigs.count(indiceNeigboor[i])==0){
			string neigboor(unitigs[indiceNeigboor[i]]);
			if (neigboor.size()>lengthRequired){
				paths.push_back(path{neigboor,indiceNeigboor[i]});
			}else{
				vector<path> paths2add;
				vector<path> paths2(listPath(lengthRequired-neigboor.size()+kgraph+1, indiceNeigboor[i],usedUnitigs));
				for (size_t j(0);j<paths2.size();++j){
					string str (compaction(paths2[j].str,unitig,kgraph));
					if(!str.empty()){
						paths2add.push_back({str,paths2[j].lastUnitig});
					}else{

					}
				}
				paths.insert(paths.end(), paths2add.begin(), paths2add.end());
			}
		}
	}
	return paths;
}


//vector<path> MappingSupervisor::listPathSons(size_t lengthRequired, uNumber ind, unordered_set<uNumber>& usedUnitigs){
//	string unitig(unitigs[ind]);
//	vector<path> paths;
//	vector<uNumber> indiceNeigboor(graph[getRepresent(unitig.substr(unitig.size()-kgraph,kgraph))]);
//	usedUnitigs.insert(ind);
//	//foreach unvisited neighbor
//	for (size_t i(0); i<indiceNeigboor.size(); ++i){
//		if(usedUnitigs.count(indiceNeigboor[i])==0){
//			string neigboor(unitigs[indiceNeigboor[i]]);
//			if (neigboor.size()>lengthRequired){
//				paths.push_back(path{neigboor,indiceNeigboor[i]});
//			}else{
//				vector<path> paths2(listPathSons(lengthRequired-neigboor.size()+kgraph, indiceNeigboor[i],usedUnitigs));
//				for (size_t j(0);j<paths2.size();++j){
//					string str (compaction(paths2[j].str,unitig,kgraph));
//					if(!str.empty()){
//						paths2[j].str=str;
//					}else{
//
//					}
//				}
//				paths.insert(paths.end(), paths2.begin(), paths2.end());
//			}
//		}
//	}
//	return paths;
//}

//
//vector<path> MappingSupervisor::listPathFathers(size_t lengthRequired, uNumber ind, unordered_set<uNumber>& usedUnitigs){
//	string unitig(unitigs[ind]);
//	vector<path> paths;
//	vector<uNumber> indiceNeigboor(graph[getRepresent(unitig.substr(0,kgraph))]);
//	vector<uNumber> indiceNeigboor2(graph[getRepresent(unitig.substr(unitig.size()-kgraph,kgraph))]);
//	indiceNeigboor.insert(indiceNeigboor.end(), indiceNeigboor2.begin(), indiceNeigboor2.end());
//	usedUnitigs.insert(ind);
//	//foreach unvisited neighbor
//	for (size_t i(0); i<indiceNeigboor.size(); ++i){
//		if(usedUnitigs.count(indiceNeigboor[i])==0){
//			string neigboor(unitigs[indiceNeigboor[i]]);
//			if (neigboor.size()>lengthRequired){
//				paths.push_back(path{neigboor,indiceNeigboor[i]});
//			}else{
//				vector<path> paths2(listPathFathers(lengthRequired-neigboor.size()+kgraph, indiceNeigboor[i],usedUnitigs));
//				for (size_t j(0);j<paths2.size();++j){
//					string str (compaction(paths2[j].str,unitig,kgraph));
//					if(!str.empty()){
//						paths2[j].str=str;
//					}else{
//
//					}
//				}
//				paths.insert(paths.end(), paths2.begin(), paths2.end());
//			}
//		}
//	}
//	return paths;
//}



void MappingSupervisor::findCandidate(const string& unitig, unordered_set<minimizer>& min, unordered_map<rNumber,size_t>& Candidate, unordered_map<rNumber,unordered_set<minimizer>>& read2Min){
	if(unitigs.size()<minSizeUnitigs){
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
				genomicKmers=allKmerSet(k2,homocompression( unitig));
			}
			if(jaccard(k2,homocompression(region),genomicKmers)>minJacc){
				readMapped++;
				return ((int)j);
			}
		}
	}
	return -1;
}



bool MappingSupervisor::alignOnPath(const path& path, const string& read, size_t position, unordered_set<uNumber>& usedUnitigsInitial){

	if(read.empty()){return false;}
	if(read.size()<=position+offset){return true;}
	int start ((int)position-(int)path.str.size()*1);

	string region(read.substr(max(0,start),path.str.size()*2));
	if(region.size()<offset){cout<<path.str.size()<<endl; return true;}

	unordered_set<minimizer> genomicKmers=allKmerSet(k2,homocompression(path.str));
//	unordered_set<minimizer> genomicKmers=allKmerSet(k2,path.str);

	if(jaccard(k2,homocompression(region),genomicKmers)>minJacc){
//	if(jaccard(k2,region,genomicKmers)>minJacc){

		auto list(listPath(offset, path.lastUnitig, usedUnitigsInitial));
		for(size_t i(0);i<list.size();++i){
			auto usedUnitigs(usedUnitigsInitial);
			if(alignOnPath(list[i], read, position+path.str.size(), usedUnitigs)){
//				outFile<<"read : "<<read<<" path:  "<<path.str<<endl;
				return true;
			}
		}
	}
	return false;
}




void MappingSupervisor::MapPart(size_t L, size_t R){
	//foreach unitig (sort of)
	for (size_t i(L); i<R; ++i){
		if(unitigsPreMapped%100==0){
			cout<<unitigsPreMapped++<<endl;
//			cout<<aligneOnPathSucess<<endl;
		}

		string unitig=unitigs[i];
		if(unitig.size()>minSizeUnitigs){
			bigUnitig++;
			unordered_set<uNumber> usedUnitigsShared;
			vector<path> list(listPath(offset, (uNumber)i, usedUnitigsShared));
			unordered_set<minimizer> min;
			bool done(false);
			unordered_map<rNumber,size_t> Candidate;Candidate.set_empty_key(-1);
			unordered_map<rNumber,unordered_set<minimizer>> read2min;read2min.set_empty_key(-1);
			unordered_set<minimizer> genomicKmers;
			findCandidate(unitig,min,Candidate,read2min);
			//foreach reads that could map on the unitig (prepremapped)
			for(auto it=Candidate.begin(); it!=Candidate.end(); ++it){
				if(it->second>=multi){
					int position(isCandidateCorrect(unitig,it->first,read2min,genomicKmers));
					if(position!=-1){
						//the read is mapped on the unitig (pre-mapped)
						if(!done){
							unitigsPreMapped++;
							done=true;
						}
						bool mappedOnGraph(false);
						//For each possible path
						for(size_t ii(0); ii<list.size() and !mappedOnGraph; ++ii){
							unordered_set<uNumber> usedUnitigs(usedUnitigsShared);
							path path(list[ii]);
							string read(reads[it->first]);
							if(alignOnPath(path, read.substr(0,position), 0, usedUnitigs)){
								if(read.size()<position+unitig.size()){
									mappedOnGraph=true;
								}else{
									if(alignOnPath(path, read.substr(position), 0, usedUnitigs)){
										mappedOnGraph=true;
									}
								}
								if(mappedOnGraph){
									++aligneOnPathSucess;
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






void MappingSupervisor::MapAll(){
	size_t nbThreads(4);
	vector<thread> threads;
	vector<size_t> limits = bounds(nbThreads, unitigs.size());

	for(size_t i(0); i<nbThreads; ++i) {
		threads.push_back(thread(&MappingSupervisor::MapPart, this, limits[i], limits[i+1]));
	}

	for(auto &t : threads){t.join();}
	cout<<"Unitigs mapped "<<unitigsPreMapped<<" Percent unitigs mapped : "<<(100*unitigsPreMapped)/(bigUnitig)<<endl;
	cout<<"Read pre-mapped : "<<readMapped<<" Percent read pre-mapped : "<<(100*readMapped)/(reads.size()+1)<<endl;
	cout<<"Read mapped on graph: "<<aligneOnPathSucess<<" Percent read mapped  on graph: "<<(100*aligneOnPathSucess)/(reads.size()+1)<<endl;
}





