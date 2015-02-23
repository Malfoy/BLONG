//
//  MappingSupervisor.cpp
//  PBMOG
//
//  Created by malfoy on 09/01/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#include "MappingSupervisor.h"
#include "Utils.h"




//vector<path> MappingSupervisor::listPathSons(size_t lengthRequired, const string& substr, char depth){
//	if(depth>depthMax){
//		return {};
//	}
//	cout<<"go"<<endl;
////	string unitig(unitigs[uNum]);
//	cout<<"sub :"<<substr<<endl;
//	vector<path> paths;
//	vector<uNumber> Sons(G.getBegin(substr));
//	//foreach son
//	for (size_t i(0); i<Sons.size(); ++i){
//		string son(unitigs[Sons[i]]);
//		cout<<"son "<<son<<endl;
//		if (son.size()>lengthRequired){
//			paths.push_back(path{son,Sons[i],{Sons[i]}});
//		}else{
//			vector<uNumber> grandSons(G.getBegin(son.substr(son.size()-kgraph,kgraph)));
//			//foreach grandsons
//			for (size_t j(0); j<grandSons.size(); ++j){
//				string songrandson (compactionEnd(son,unitigs[grandSons[j]],kgraph));
//				cout<<"grandson"<<unitigs[grandSons[j]]<<" "<<grandSons[j]<<endl;
//				cout<<"songrandson "<<songrandson<<endl;
//				if(!songrandson.empty()){
//					if(lengthRequired>songrandson.size()){
//						vector<path> recurPaths(listPathSons(lengthRequired-songrandson.size()+kgraph, songrandson.substr(songrandson.size()-kgraph, kgraph), depth+1));
//						for (size_t ii(0);ii<recurPaths.size();++ii){
//							string sonGrandSonAndCo(compactionEnd(songrandson, recurPaths[ii].str, kgraph));
//							cout<<"recurpath "<<recurPaths[ii].str<<endl;
//
//							if(!sonGrandSonAndCo.empty()){
//								recurPaths[ii].numbers.push_back(Sons[i]);
//								recurPaths[ii].numbers.push_back(grandSons[j]);
//
//								paths.push_back({sonGrandSonAndCo, recurPaths[ii].lastUnitig, recurPaths[ii].numbers});
//							}else{
//								cout<<"wtfson"<<endl;
//								cout<<songrandson<<" "<<recurPaths[ii].str<<endl;
//
//								cin.get();
//							}
//						}
//					}else{
//						paths.push_back({songrandson, grandSons[j], {Sons[i],grandSons[j]}});
//					}
//				}else{
//					cout<<"wtf"<<endl;
//					cin.get();
//
//				}
//			}
//		}
//	}
//	cout<<"end"<<endl;
//	return paths;
//}

//vector<path> MappingSupervisor::listPathSons(size_t lengthRequired,uint32_t ind, char depth){
//	if(depth>depthMax){return {};}
//
//	vector<path> result;
//	const string unitig(unitigs[ind]);
//	string beg(unitig.substr(0,kgraph)),end(unitig.substr(unitig.size()-kgraph,kgraph));
//	string begRC(reversecomplement(beg)),endRC(reversecomplement(end));
//
//	vector<uNumber> sons,fathers;
//	if(beg<begRC){
//		sons=(G.getEnd(beg));
//	}else{
//		sons=(G.getBegin(begRC));
//	}
//	if(end<endRC){
//		fathers=(G.getBegin(end));
//	}else{
//		fathers=(G.getEnd(endRC));
//	}
//
//	for(size_t i(0);i<sons.size();++i){
//
//	}
//}



vector<path> MappingSupervisor::listPathSons(size_t lengthRequired, const string& substr, char depth){
//	cout<<substr<<endl;
//	if(depth>depthMax){return {};}
	vector<path> paths;
	vector<uNumber> Sons(G.getBegin(substr));
//	foreach son
	for (size_t i(0); i<Sons.size(); ++i){
		string son(unitigs[Sons[i]]);
//		cout<<"son :"<<son<<endl;
		string newson(compactionEnd(substr, son, kgraph));
//		if(newson.substr(newson.size()-kgraph)==substr){
//			cout<<"omg"<<endl;cin.get();
//		}
		if(newson.empty()){
			cout<<"wtfnewson"<<endl;
			cin.get();
		}
		if (son.size()>=lengthRequired){
			paths.push_back(path{son,{Sons[i]}});
		}else{
			vector<path> recurPaths(listPathSons(lengthRequired-son.size()+kgraph, newson.substr(newson.size()-kgraph), depth+1));
			for (size_t j(0); j<recurPaths.size(); ++j){
//				cout<<"recur"<<endl;
				//should be compactionEnd
				string sonAndCo(compactionEnd(newson,recurPaths[j].str, kgraph));
				if(sonAndCo.size()>=lengthRequired){
//					cout<<"sonAndCo : "<<sonAndCo<<endl;
					recurPaths[j].numbers.push_back(Sons[i]);
					paths.push_back(path{sonAndCo, recurPaths[j].numbers});
				}else{
					cout<<lengthRequired<<endl;
					cout<<sonAndCo.size()<<endl;
					cout<<"wtfson"<<endl;
					cout<<son<<" "<<recurPaths[j].str<<endl;
					cout<<son.substr(0,kgraph)<<" "<<getRepresent(son.substr(son.size()-kgraph))<<endl;;
					cout<<recurPaths[j].str.substr(0,kgraph)<<" "<<getRepresent(recurPaths[j].str.substr(recurPaths[j].str.size()-kgraph))<<endl;;
					cin.get();
				}
			}
		}
	}
return paths;
}



vector<path> MappingSupervisor::listPathFathers(size_t lengthRequired, const string& substr, char depth){
//	if(depth>depthMax){return {};}
	vector<path> paths;
	vector<uNumber> Fathers(G.getEnd(substr));
	//foreach son
	for (size_t i(0); i<Fathers.size(); ++i){
		string father(unitigs[Fathers[i]]);
		string newfather(compactionBegin(substr, father, kgraph));
		if(newfather.empty()){
			cout<<"wtfnewfather"<<endl;
			cin.get();
		}

		if (father.size()>=lengthRequired){
			paths.push_back(path{father,{Fathers[i]}});
		}else{
			vector<path> pathsRecur(listPathFathers(lengthRequired-father.size()+kgraph,newfather.substr(0,kgraph), depth+1));
			for (size_t j(0); j<pathsRecur.size(); ++j){
				//should be compactionbegin
				string fatherandco(compactionBegin(newfather,pathsRecur[j].str,kgraph));
				if(fatherandco.size()>=lengthRequired){
					pathsRecur[j].numbers.push_back(Fathers[i]);
					paths.push_back({fatherandco, pathsRecur[j].numbers});
				}else{
					cout<<father<<" "<<pathsRecur[j].str<<endl;
					cout<<"wtffather"<<endl;
					cin.get();
				}
			}
		}
	}
	return paths;
}



void MappingSupervisor::findCandidate(const string& unitig, unordered_set<minimizer>& min, unordered_map<rNumber,size_t>& Candidate, unordered_map<rNumber,unordered_set<minimizer>>& read2Min){
	if(unitigs.size()<H){
		for(size_t j(0);j+k<unitig.size();++j){
			minimizer seq=seq2int(unitig.substr(j,k));
			if(min.count(seq)==0){
				min.insert(seq);
				if(min2Reads.count(seq)!=0){
					vector<rNumber> myset=min2Reads[seq];
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
					vector<rNumber> myset=min2Reads[seq];
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
			if(jaccard(k2,region,genomicKmers)>=minJacc){
				readMapped++;
				return ((int)j);
			}
		}
	}
	return -1;
}




bool MappingSupervisor::alignOnPathSons(const path& path, const string& read, size_t position,vector<uNumber>& numbers){
	if(read.empty()){return false;}
	if(read.size()<=position+offset){return true;}
	int start ((int)position-(int)path.str.size()*1);

	string region(read.substr(max(0,start),path.str.size()*2));
	if(region.size()<offset){return true;}


//		unordered_set<minimizer> genomicKmers=allKmerSet(k2,homocompression(path.str));
	unordered_set<minimizer> genomicKmers=allKmerSet(k2,path.str);

//		if(jaccard(k2,homocompression(region),genomicKmers)>=minJacc){
	if(jaccard(k2,read,genomicKmers)>=minJacc){
		regionmapped++;
		if(read.size()<=position+path.str.size()+offset){
			return true;
		}else{
			size_t size(numbers.size());
			numbers.insert(numbers.end(),path.numbers.begin(),path.numbers.end());
			auto list(listPathSons(offset, path.str.substr(path.str.size()-kgraph,kgraph),0));
			for(size_t i(0);i<list.size();++i){
				if(alignOnPathSons(list[i], read, position+path.str.size(),numbers)){
					return true;
				}
			}
			numbers.resize(size);
		}
	}
	return false;
}


bool MappingSupervisor::alignOnPathFathers(const path& path, const string& read, size_t position,vector<uNumber>& numbers ){
//	cout<<"father : "<<path.str<<endl;
	if(read.empty()){return false;}
	if(read.size()<=position+offset){return true;}
	int start ((int)position-(int)path.str.size()*1);

	string region(read.substr(max(0,start),path.str.size()*2));
	if(region.size()<offset){return true;}

//		unordered_set<minimizer> genomicKmers=allKmerSet(k2,homocompression(path.str));
	unordered_set<minimizer> genomicKmers=allKmerSet(k2,path.str);

//		if(jaccard(k2,homocompression(region),genomicKmers)>=minJacc){
	if(jaccard(k2,read,genomicKmers)>=minJacc){
		regionmapped++;
		if(read.size()<=position+path.str.size()+offset){
			return true;
		}else{
			size_t size(numbers.size());
			numbers.insert(numbers.end(),path.numbers.begin(),path.numbers.end());
			auto list(listPathFathers(offset, path.str.substr(0,kgraph),0));
			for(size_t i(0);i<list.size();++i){
				if(alignOnPathFathers(list[i], read, position+path.str.size(),numbers)){
					return true;
				}
			}
			numbers.resize(size);
		}
	}
	return false;
}




void MappingSupervisor::MapPart(size_t L, size_t R){
	//foreach unitig (sort of)
	for (size_t i(L); i<R; ++i){
		string unitig=unitigs[i];
//		cout<<"unitig"<<unitig<<endl;
		if(unitig.size()>minSizeUnitigs){
			if(bigUnitig++%1000==0){
				cout<<bigUnitig<<endl;
				cout<<"Unitigs mapped "<<unitigsPreMapped<<" Percent unitigs mapped : "<<(100*unitigsPreMapped)/(bigUnitig+1)<<endl;
				cout<<"Read pre-mapped : "<<readMapped<<" Percent read pre-mapped : "<<(100*readMapped)/(reads.size()+1)<<endl;
				cout<<"Read mapped on graph: "<<aligneOnPathSucess<<" Percent read mapped  on graph: "<<(100*aligneOnPathSucess)/(reads.size()+1)<<endl;
				cout<<island<<" islands..."<<endl;
				cout<<regionmapped<<" region mapped"<<endl;
				cout<<endl;
			}
			vector<path> list(listPathSons(offset, unitigs[i].substr(unitigs[i].size()-kgraph,kgraph),0));
//			cout<<"nb sons"<<list.size()<<endl;
			vector<path> list2(listPathFathers(offset,unitigs[i].substr(0,kgraph),0));
			list.insert(list.end(),list2.begin(),list2.end());
			if(list.size()==0){
				island++;
//				cin.get();
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
						string read(reads[it->first]);
						//here we want to find a position
						vector<uNumber> ref;
						for(size_t ii(0); ii<list.size() and !mappedOnGraph; ++ii){
							path path(list[ii]);
							path.str=compaction(path.str, unitig, kgraph);
							if(!path.str.empty()){
//								cout<<"go"<<endl;
//								cout<<path.str<<endl;
								if(alignOnPathSons(path, read, 0,ref) or alignOnPathFathers(path, read, 0,ref)){
									mappedOnGraph=true;
								}
//								cout<<"end"<<endl;
								if(mappedOnGraph){
									++aligneOnPathSucess;
//									cout<<"go1"<<endl;
//									for_each(ref.begin(), ref.end(), [](uNumber u){cout<<u<<" ";});
//									cout<<" end"<<endl;cin.get();
									mutexEraseReads.lock();
									reads[it->first].clear();
									mutexEraseReads.unlock();
								}
							}
						}
					}
				}
			}else{
				//foreach reads that could map on the unitig (prepremapped)
				for(auto it=Candidate.begin(); it!=Candidate.end(); ++it){
//					cout<<"go3"<<endl;
					if(it->second>=multi){
//						cout<<"go4"<<endl;
						int position(isCandidateCorrect(unitig,it->first,read2min,genomicKmers));
						if(position!=-1){
							//the read is mapped on the unitig (pre-mapped)
							if(!done){
								unitigsPreMapped++;
								done=true;
							}
							bool mappedOnGraph(false);
							string read(reads[it->first]);
							vector<uNumber> ref;
							//For each possible path
							for(size_t j(0); j<list.size() and !mappedOnGraph; ++j){
								path path(list[j]);
//								cout<<"go2"<<endl;
								if(alignOnPathSons(path, read.substr(0,position), 0,ref) or alignOnPathFathers(path, read.substr(0,position), 0,ref) ){
//									cout<<"suceed"<<position<<endl;
									if(read.size()<position+unitig.size()){
										mappedOnGraph=true;
									}else{
										if(alignOnPathSons(path, read.substr(position), 0,ref) or alignOnPathFathers(path, read.substr(position), 0,ref)){
											mappedOnGraph=true;
										}
									}
									if(mappedOnGraph){
										++aligneOnPathSucess;
//										cout<<"go2"<<endl;
//										for_each(ref.begin(), ref.end(), [](uNumber u){cout<<u<<" ";});
//										cout<<" end"<<endl;cin.get();

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
	cout<<island<<" islands..."<<endl;
}





