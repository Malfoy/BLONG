//
//  MappingSupervisor.cpp
//  PBMOG
//
//  Created by malfoy on 09/01/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#include "MappingSupervisor.h"
#include "Utils.h"





vector<path> MappingSupervisor::listPathSons(size_t lengthRequired, const string& substr, char depth){
	//	cout<<substr<<endl;
	if(depth>depthMax){return {};}
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
	if(depth>depthMax){return {};}
	vector<path> paths;
	vector<uNumber> Fathers(G.getEnd(substr));
	//foreach son
	//	cout<<Fathers.size()<<endl;

	for (size_t i(0); i<Fathers.size(); ++i){
		string father(unitigs[Fathers[i]]);
		string newfather(compactionBegin(substr, father, kgraph));
		//		cout<<newfather<<endl;
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
			if(min.unordered_set::count(seq)==0){
				min.insert(seq);
				if(min2Reads.unordered_map::count(seq)!=0){
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
			if(min.unordered_set::count(seq)==0){
				min.insert(seq);
				if(min2Reads.unordered_map::count(seq)!=0){
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
		if(setMin.unordered_set::count(seq)!=0){
			string region;

			if(j>1*unitig.size()){
				region=read.substr(j-1*unitig.size(),2*unitig.size());
			}else{
				region=read.substr(0,2*unitig.size());
			}
			if(genomicKmers.size()==0){
				genomicKmers=allKmerSetStranded(k2,unitig);
			}
			//if(jaccard(k2,homocompression(region),genomicKmers)>minJacc){
			if(jaccardStranded(k2,region,genomicKmers)>=minJacc){
				readMapped++;
				return ((int)j);
			}
		}
	}
	return -1;
}

int MappingSupervisor::isCandidateCorrectReverse(const string& unitig, rNumber readNumber, unordered_map<rNumber,unordered_set<minimizer>>& read2min, unordered_set<minimizer>& genomicKmers){
	bool goodreadb=false;
	unordered_set<minimizer> setMin(read2min[readNumber]);
	string read=reversecomplement(reads[readNumber]);
	if(read.empty()){return -1;}
	for(size_t j(0);j+k<read.size() and !goodreadb;++j){
		minimizer seq(seq2int(read.substr(j,k)));
		if(setMin.unordered_set::count(seq)!=0){
			string region;

			if(j>1*unitig.size()){
				region=read.substr(j-1*unitig.size(),2*unitig.size());
			}else{
				region=read.substr(0,2*unitig.size());
			}
			if(genomicKmers.size()==0){
				genomicKmers=allKmerSetStranded(k2,unitig);
			}
			//if(jaccard(k2,homocompression(region),genomicKmers)>minJacc){
			if(jaccardStranded(k2,region,genomicKmers)>=minJacc){
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
	if(jaccard(k2,region,genomicKmers)>=minJacc){
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
	int start ((int)read.size()-(int)position-(int)path.str.size()*1);

	string region(read.substr(max(0,start),path.str.size()*2));
	if(region.size()<offset){return true;}

	//		unordered_set<minimizer> genomicKmers=allKmerSet(k2,homocompression(path.str));
	unordered_set<minimizer> genomicKmers=allKmerSet(k2,path.str);

	//		if(jaccard(k2,homocompression(region),genomicKmers)>=minJacc){
	if(jaccard(k2,region,genomicKmers)>=minJacc){
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
		//		cout<<"unitig :"<<unitig<<endl;
		if(unitig.size()>=minSizeUnitigs){
			if(bigUnitig++%1000==0){
				cout<<bigUnitig++<<endl;
				cout<<"Unitigs mapped "<<unitigsPreMapped<<" Percent unitigs mapped : "<<(100*unitigsPreMapped)/(bigUnitig+1)<<endl;
				cout<<"Read pre-mapped : "<<readMapped<<" Percent read pre-mapped : "<<(100*readMapped)/(reads.size()+1)<<endl;
				cout<<"Read mapped on graph: "<<aligneOnPathSucess<<" Percent read mapped  on graph: "<<(100*aligneOnPathSucess)/(reads.size()+1)<<endl;
				cout<<island<<" islands..."<<endl;
				cout<<regionmapped<<" region mapped"<<endl;
				cout<<leftmap<<" "<< rightmap<<endl;
				cout<<endl;
			}
			vector<path> ListSons(listPathSons(offset, unitigs[i].substr(unitigs[i].size()-kgraph,kgraph),0));
			vector<path> ListFathers(listPathFathers(offset, unitigs[i].substr(0,kgraph),0));
			if( ListFathers.size()==0 and ListSons.size()==0){
				//				cin.get();
				island++;
				continue;
			}
			unordered_set<minimizer> min;
			bool done(false);
			unordered_map<rNumber,size_t> Candidate;Candidate.set_empty_key(-1);
			unordered_map<rNumber,unordered_set<minimizer>> read2min;read2min.set_empty_key(-1);
			unordered_set<minimizer> genomicKmers,genomicKmersRC;
			findCandidate(unitig,min,Candidate,read2min);
			//			cout<<1<<endl;
			//			if(unitig.size()<offset){
			if(/* DISABLES CODE */ (false)){
				for(auto it=Candidate.begin(); it!=Candidate.end(); ++it){
					if(it->second>=multi){
						bool mappedOnGraph(false);
						string read(reads[it->first]);
						//here we want to find a position
						vector<uNumber> ref;
						for(size_t ii(0); ii<ListSons.size() and !mappedOnGraph; ++ii){
							path path(ListSons[ii]);
							path.str=compaction(path.str, unitig, kgraph);
							if(!path.str.empty()){
								//								cout<<"go"<<endl;
								//								cout<<path.str<<endl;
								if(alignOnPathSons(path, read, 0,ref) or alignOnPathFathers(path, read, 0,ref)){
									mappedOnGraph=true;
								}

								if(mappedOnGraph){
									++aligneOnPathSucess;
									cout<<unitig<<endl;
									//									cout<<"go1"<<endl;
									//									for_each(ref.begin(), ref.end(), [](uNumber u){cout<<u<<" ";});
									//									cout<<" end"<<endl;cin.get();
									//									mutexEraseReads.lock();
									//									reads[it->first].clear();
									//									mutexEraseReads.unlock();
								}
							}
						}
					}
				}
			}else{

				//foreach reads that could map on the unitig (prepremapped)
				for(auto it=Candidate.begin(); it!=Candidate.end(); ++it){
					if(it->second>=multi){
						int position;
						int position1(isCandidateCorrect(unitig,it->first,read2min,genomicKmers));
						int position2(isCandidateCorrectReverse(unitig,it->first,read2min,genomicKmersRC));
						string read;
						if(position1!=-1){
							position=position1;
							read=(reads[it->first]);
						}else{
							if (position2!=-1) {
								position=position2;
								read=reversecomplement(reads[it->first]);
							}else{
								continue;
							}
						}
						//the read is mapped on the unitig (pre-mapped)
						if(!done){
							unitigsPreMapped++;
							done=true;
						}
						bool mappedLeft(false),mappedRight(false);

						vector<uNumber> ref;
						//For each possible path
						if(position<(int)offset){
							mappedLeft=true;
							leftmap++;
						}else{
							for(size_t j(0); j<ListFathers.size(); ++j){
								path pathFather(ListFathers[j]);

								if(alignOnPathFathers(pathFather, read.substr(0,position), 0,ref)){
									mappedLeft=true;
									leftmap++;
									break;
								}
							}
						}
						if(read.size()<position+unitig.size()+offset){
							mappedRight=true;
							rightmap++;
						}else{
							for(size_t j(0); j<ListSons.size(); ++j){
								path pathSon(ListSons[j]);
								if(alignOnPathSons(pathSon, read.substr(position+unitig.size()), 0,ref)){
									mappedRight=true;
									rightmap++;
									break;
								}
							}
						}

						if(mappedRight and mappedLeft){
							++aligneOnPathSucess;
//							cout<<"go2"<<endl;
//							for_each(ref.begin(), ref.end(), [this](uNumber u){cout<<unitigs[u]<<" ";});
//							cout<<" end"<<endl;cin.get();

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









void MappingSupervisor::MapAll(){
	size_t nbThreads(4);
	vector<thread> threads;
	vector<size_t> limits = bounds(nbThreads, unitigs.size());
	cout<<unitigs.size()<<endl;

	for(size_t i(0); i<nbThreads; ++i) {
		threads.push_back(thread(&MappingSupervisor::MapPart, this, limits[i], limits[i+1]));
	}

	for(auto &t : threads){t.join();}
	cout<<"Unitigs mapped "<<unitigsPreMapped<<" Percent unitigs mapped : "<<(100*unitigsPreMapped)/(bigUnitig+1)<<endl;
	cout<<"Read pre-mapped : "<<readMapped<<" Percent read pre-mapped : "<<(100*readMapped)/(reads.size()+1)<<endl;
	cout<<"Read mapped on graph: "<<aligneOnPathSucess<<" Percent read mapped  on graph: "<<(100*aligneOnPathSucess)/(reads.size()+1)<<endl;
	cout<<island<<" islands..."<<endl;
}





