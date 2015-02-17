//
//  MappingSupervisor.cpp
//  PBMOG
//
//  Created by malfoy on 09/01/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#include "MappingSupervisor.h"
#include "Utils.h"




//
//vector<path> MappingSupervisor::listPath(size_t lengthRequired, uNumber ind, set<uNumber> usedUnitigs){
//	string unitig(unitigs[ind]);
//	vector<path> paths;
//	vector<uNumber> indiceNeigboor(G.getEnd(unitig.substr(0,kgraph)));
//	vector<uNumber> indiceNeigboor2(G.getBegin(unitig.substr(unitig.size()-kgraph,kgraph)));
////	cout<<indiceNeigboor.size()<<" "<<indiceNeigboor2.size()<<endl;
////	cout<<getRepresent(unitig.substr(0,kgraph))<<endl;
////	for(size_t i(0);i<indiceNeigboor.size();++i){
////		cout<<unitigs[indiceNeigboor[i]]<<endl;
////
////	}
////	cin.get();
//	indiceNeigboor.insert(indiceNeigboor.end(), indiceNeigboor2.begin(), indiceNeigboor2.end());
//	usedUnitigs.insert(ind);
//	//foreach unvisited neighbor
//	for (size_t i(0); i<indiceNeigboor.size(); ++i){
//		if(usedUnitigs.count(indiceNeigboor[i])==0){
//			string neigboor(unitigs[indiceNeigboor[i]]);
//			if (neigboor.size()>lengthRequired){
//				string str(compaction(unitig,neigboor,kgraph));
//				if(str!=""){
//					paths.push_back(path{str,indiceNeigboor[i]});
//				}else{
//					cout<<"nadine"<<endl;
//					cin.get();
//				}
//			}else{
//				vector<path> paths2add;
//				vector<path> paths2(listPath(lengthRequired-neigboor.size()+kgraph+1, indiceNeigboor[i],usedUnitigs));
//				for (size_t j(0);j<paths2.size();++j){
//					string str (compaction(paths2[j].str,neigboor,kgraph));
//					if(!str.empty()){
//						paths2add.push_back({str,paths2[j].lastUnitig});
////						cout<<"it works"<<endl;
//					}else{
////						cout<<""wut"<<endl;
//					}
//				}
//				paths.insert(paths.end(), paths2add.begin(), paths2add.end());
//			}
//		}
//	}
//	return paths;
//}


vector<path> MappingSupervisor::listPathSons(size_t lengthRequired, uNumber ind, char depth){
	if(depth>depthMax){
//		cout<<"max"<<endl;
		return {};
	}
	string unitig(unitigs[ind]);
	vector<path> paths;
	vector<uNumber> Sons(G.getBegin(unitig.substr(unitig.size()-kgraph,kgraph)));
	//foreach son
	for (size_t i(0); i<Sons.size(); ++i){
		string son(unitigs[Sons[i]]);
		if (son.size()>lengthRequired){
			paths.push_back(path{son,Sons[i]});
		}else{
			vector<uNumber> grandSons(G.getBegin(son.substr(son.size()-kgraph,kgraph)));
			//foreach grandsons
			for (size_t j(0); j<grandSons.size(); ++j){
				string str (compaction(unitigs[grandSons[j]],son,kgraph));
				if(!str.empty()){
					if(lengthRequired-str.size()>0){
						vector<path> paths2(listPathSons(lengthRequired-str.size(), grandSons[j],depth+1));
						for (size_t ii(0);ii<paths2.size();++ii){
							string str2(compaction(str,paths2[ii].str,kgraph));
							if(!str.empty()){
								paths.push_back({str2,paths2[ii].lastUnitig});
							}else{
								cout<<"wtf"<<endl;
							}
						}
					}
				}else{
					cout<<"wtf"<<endl;
				}
			}
		}
	}
	return paths;
}



vector<path> MappingSupervisor::listPathFathers(size_t lengthRequired, uNumber ind, char depth){
	if(depth>depthMax){
//		cout<<"max"<<endl;
		return {};
	}
	string unitig(unitigs[ind]);
	vector<path> paths;
	vector<uNumber> Sons(G.getEnd(unitig.substr(0,kgraph)));
	//foreach son
	for (size_t i(0); i<Sons.size(); ++i){
		string son(unitigs[Sons[i]]);
		if (son.size()>lengthRequired){
			paths.push_back(path{son,Sons[i]});
		}else{
			vector<uNumber> grandSons(G.getEnd(son.substr(0,kgraph)));
			//foreach grandsons
			for (size_t j(0); j<grandSons.size(); ++j){
				string str (compaction(unitigs[grandSons[j]],son,kgraph));
				if(!str.empty()){
					if(lengthRequired-str.size()>0){
						vector<path> paths2(listPathSons(lengthRequired-str.size(), grandSons[j],depth+1));
						for (size_t ii(0);ii<paths2.size();++ii){
							string str2(compaction(str,paths2[ii].str,kgraph));
							if(!str.empty()){
								paths.push_back({str2,paths2[ii].lastUnitig});
							}else{
								cout<<"wtf"<<endl;
							}
						}
					}
				}else{
					cout<<"wtf"<<endl;
				}
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
		for(size_t j(0);j<sketch.size();++j){
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
			//if(jaccard(k2,homocompression(region),genomicKmers)>minJacc){
			if(jaccard(k2,region,genomicKmers)>minJacc){
				readMapped++;
				return ((int)j);
			}
		}
	}
	return -1;
}



bool MappingSupervisor::alignOnPath(const path& path, const string& read, size_t position, set<uNumber>& usedUnitigsInitial){

	if(read.empty()){return false;}
	if(read.size()<=position+offset){return true;}
	int start ((int)position-(int)path.str.size()*1);

	string region(read.substr(max(0,start),path.str.size()*2));
	if(region.size()<offset){
//		cout<<path.str.size()<<endl;
		return true;}

//	unordered_set<minimizer> genomicKmers=allKmerSet(k2,homocompression(path.str));
	unordered_set<minimizer> genomicKmers=allKmerSet(k2,path.str);

//	if(jaccard(k2,homocompression(region),genomicKmers)>minJacc){
	if(jaccard(k2,region,genomicKmers)>minJacc){

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


bool MappingSupervisor::alignOnPathSons(const path& path, const string& read, size_t position){

	if(read.empty()){return false;}
	if(read.size()<=position+offset){return true;}
	int start ((int)position-(int)path.str.size()*1);

	string region(read.substr(max(0,start),path.str.size()*2));
	if(region.size()<offset){
//		cout<<path.str.size()<<endl;
		return true;}

	//	unordered_set<minimizer> genomicKmers=allKmerSet(k2,homocompression(path.str));
	unordered_set<minimizer> genomicKmers=allKmerSet(k2,path.str);

	//	if(jaccard(k2,homocompression(region),genomicKmers)>minJacc){
	if(jaccard(k2,read,genomicKmers)>minJacc){

		auto list(listPathSons(offset, path.lastUnitig,0));
		for(size_t i(0);i<list.size();++i){
			if(alignOnPathSons(list[i], read, position+path.str.size())){
				//				outFile<<"read : "<<read<<" path:  "<<path.str<<endl;
				return true;
			}
		}
	}
	return false;
}


bool MappingSupervisor::alignOnPathFathers(const path& path, const string& read, size_t position){

	if(read.empty()){return false;}
	if(read.size()<=position+offset){return true;}
	int start ((int)position-(int)path.str.size()*1);

	string region(read.substr(max(0,start),path.str.size()*2));
	if(region.size()<offset){
//		cout<<path.str.size()<<endl;
		return true;}

	//	unordered_set<minimizer> genomicKmers=allKmerSet(k2,homocompression(path.str));
	unordered_set<minimizer> genomicKmers=allKmerSet(k2,path.str);

	//	if(jaccard(k2,homocompression(region),genomicKmers)>minJacc){
	if(jaccard(k2,read,genomicKmers)>minJacc){

		auto list(listPathFathers(offset, path.lastUnitig,0));
		for(size_t i(0);i<list.size();++i){
			if(alignOnPathFathers(list[i], read, position+path.str.size())){
				//outFile<<"read : "<<read<<" path:  "<<path.str<<endl;
				return true;
			}
		}
	}
	return false;
}




void MappingSupervisor::MapPart(size_t L, size_t R){
	//foreach unitig (sort of)
	for (size_t i(L); i<R; ++i){
		string unitig=unitigs[i];
		if(unitig.size()>minSizeUnitigs){
			if(bigUnitig++%100==0){
				cout<<bigUnitig<<endl;
				cout<<"Unitigs mapped "<<unitigsPreMapped<<" Percent unitigs mapped : "<<(100*unitigsPreMapped)/(bigUnitig+1)<<endl;
				cout<<"Read pre-mapped : "<<readMapped<<" Percent read pre-mapped : "<<(100*readMapped)/(reads.size()+1)<<endl;
				cout<<"Read mapped on graph: "<<aligneOnPathSucess<<" Percent read mapped  on graph: "<<(100*aligneOnPathSucess)/(reads.size()+1)<<endl;
				cout<<island<<" islands..."<<endl;
				cout<<endl;
			}
			set<uNumber> usedUnitigsShared;
			vector<path> list(listPathSons(offset, (uNumber)i,0));
			vector<path> list2(listPathFathers(offset, (uNumber)i,0));
			list.insert(list.end(),list2.begin(),list2.end());
			if(list.size()==0){
				island++;
				continue;
			}
			unordered_set<minimizer> min;
			bool done(false);
			unordered_map<rNumber,size_t> Candidate;Candidate.set_empty_key(-1);
			unordered_map<rNumber,unordered_set<minimizer>> read2min;read2min.set_empty_key(-1);
			unordered_set<minimizer> genomicKmers;
			findCandidate(unitig,min,Candidate,read2min);
			if(unitig.size()<offset){
				for(auto it=Candidate.begin(); it!=Candidate.end(); ++it){
					if(it->second>=multi){
						bool mappedOnGraph(false);
						for(size_t ii(0); ii<list.size() and !mappedOnGraph; ++ii){
							path path(list[ii]);
							path.str=compaction(path.str, unitig, kgraph);
							string read(reads[it->first]);
							if(alignOnPathSons(path, read, 0) or alignOnPathFathers(path, read, 0)){
								mappedOnGraph=true;
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
			}else{
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
								path path(list[ii]);
								string read(reads[it->first]);
								if(alignOnPathSons(path, read.substr(0,position), 0) or alignOnPathFathers(path, read.substr(0,position), 0) ){
									if(read.size()<position+unitig.size()){
										mappedOnGraph=true;
									}else{
										if(alignOnPathSons(path, read.substr(position), 0) or alignOnPathFathers(path, read.substr(position), 0)){
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
}






void MappingSupervisor::MapAll(){
	size_t nbThreads(3);
	vector<thread> threads;
	vector<size_t> limits = bounds(nbThreads, unitigs.size());

	for(size_t i(0); i<nbThreads; ++i) {
		threads.push_back(thread(&MappingSupervisor::MapPart, this, limits[i], limits[i+1]));
	}

	for(auto &t : threads){t.join();}
	cout<<"Unitigs mapped "<<unitigsPreMapped<<" Percent unitigs mapped : "<<(100*unitigsPreMapped)/(bigUnitig)<<endl;
	cout<<"Read pre-mapped : "<<readMapped<<" Percent read pre-mapped : "<<(100*readMapped)/(reads.size()+1)<<endl;
	cout<<"Read mapped on graph: "<<aligneOnPathSucess<<" Percent read mapped  on graph: "<<(100*aligneOnPathSucess)/(reads.size()+1)<<endl;
	cout<<island<<" islands..."<<endl;
}





