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
			//			cout<<newson<<"newson"<<endl;
			paths.push_back(path{newson,{Sons[i]}});
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
	//	cout<<"listPathFathers"<<endl;
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
			paths.push_back(path{newfather,{Fathers[i]}});
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



bool MappingSupervisor::isCandidateCorrect(const string& unitig, rNumber readNumber, unordered_map<rNumber,unordered_set<minimizer>>& read2min, unordered_set<minimizer>& genomicKmers,int& position){
	unordered_set<minimizer> setMin(read2min[readNumber]);
	string read=reads[readNumber];
	//	cout<<unitig<<" "<<read<<endl;
	if(read.empty()){return false;}
	for(size_t j(0);j+k<=read.size();++j){
		//		cout<<read.substr(j,k)<<endl;
		minimizer seq(seq2int(read.substr(j,k)));
		int pos((int)j-(int)positionInSeq(unitig, seq, k));
		if(setMin.unordered_set::count(seq)!=0){
			string region;

			if(pos>=0){
				region=read.substr(pos,1*unitig.size());
			}else{
				region=read.substr(0,unitig.size());
			}
			if(genomicKmers.size()==0){
				genomicKmers=allKmerSetStranded(k2,unitig);
			}
			//if(jaccard(k2,homocompression(region),genomicKmers)>minJacc){
			//			cout<<region<<" "<<unitig<<endl;
			if(jaccardStranded(k2,region,genomicKmers)>=minJacc){
				//				cout<<region<<" "<<unitig<<endl;
				readMapped++;
				position=(pos);
				//				cout<<"true"<<endl;
				return true;

			}
		}
	}
	//	cout<<"false"<<endl;
	return false;
}



bool MappingSupervisor::isCandidateCorrectReverse(const string& unitig, rNumber readNumber, unordered_map<rNumber,unordered_set<minimizer>>& read2min, unordered_set<minimizer>& genomicKmers,int& position){
	unordered_set<minimizer> setMin(read2min[readNumber]);
	string read=reversecomplement(reads[readNumber]);
	//	cout<<unitig<<" "<<read<<endl;
	if(read.empty()){return false;}
	for(size_t j(0);j+k<=read.size();++j){
		minimizer seq(seq2int(read.substr(j,k)));

		if(setMin.unordered_set::count(seq)!=0){
			//			cout<<j<<" "<<positionInSeq(unitig, seq, k)<<endl;
			int pos((int)j-(int)positionInSeq(unitig, seq, k));
			string region;

			if(pos>=0){
				region=read.substr(pos,1*unitig.size());
			}else{
				region=read.substr(0,unitig.size());
			}
			if(genomicKmers.size()==0){
				genomicKmers=allKmerSetStranded(k2,unitig);
			}
			//if(jaccard(k2,homocompression(region),genomicKmers)>minJacc){
			//			cout<<region<<" "<<unitig<<endl;
			if(jaccardStranded(k2,region,genomicKmers)>=minJacc){
				//				cout<<region<<" "<<unitig<<endl;
				readMapped++;
				position=(pos);
				return true;

			}
		}
	}
	return false;
}




bool MappingSupervisor::alignOnPathSons2(const path& path, const string& read, size_t position,vector<uNumber>& numbers){
	//	if(read.size()<=position+offset){return true;}
	int start ((int)position-(int)path.str.size()*1);

	string region(read.substr(max(0,start),path.str.size()*2));
	//	if(region.size()<offset){return true;}


	//		unordered_set<minimizer> genomicKmers=allKmerSet(k2,homocompression(path.str));
	unordered_set<minimizer> genomicKmers=allKmerSet(k2,path.str);
	//	cout<<region<<" "<<path.str<<endl;
	//		if(jaccard(k2,homocompression(region),genomicKmers)>=minJacc){
	if(jaccard(k2,region,genomicKmers)>=minJacc){
		//		cout<<"good"<<endl;
		regionmapped++;
		if(read.size()<position+path.str.size()-kgraph+offset){
			return true;
		}else{
			size_t size(numbers.size());
			numbers.insert(numbers.end(),path.numbers.begin(),path.numbers.end());
			auto list(listPathSons(offset, path.str.substr(path.str.size()-kgraph,kgraph),0));
			for(size_t i(0);i<list.size();++i){
				if(alignOnPathSons(list[i], read, position+path.str.size()-kgraph,numbers)){
					return true;
				}
			}
			numbers.resize(size);
		}
	}
	return false;
}


bool MappingSupervisor::alignOnPathSons(const path& path, const string& read, size_t position,vector<uNumber>& numbers){
//	cout<<"go"<<endl;
	unordered_set<minimizer> genomicKmers=allKmerSet(k2,path.str);
	minimizer kmer;
	int start(0);
	bool found(false);
	for(uint i(0); i<path.str.size(); ++i){
		if(read.size()>position+i+k2){
			kmer=seq2int(read.substr(position+i,k2));
			if(genomicKmers.unordered_set::count(kmer)!=0){
				start=i+(int)position-(int)positionInSeq(path.str, kmer, k2);
				found=true;
				break;
			}
		}
	}
//	cout<<"found"<<endl;
	if(found){
		string region(read.substr(max(0,start),path.str.size()));
		if(jaccard(k2,region,genomicKmers)>=minJacc){
			regionmapped++;
			if(read.size()<position+path.str.size()-kgraph+offset){
				return true;
			}else{
				size_t size(numbers.size());
				numbers.insert(numbers.end(),path.numbers.begin(),path.numbers.end());
				auto list(listPathSons(offset, path.str.substr(path.str.size()-kgraph,kgraph),0));
				for(size_t i(0);i<list.size();++i){
					if(alignOnPathSons(list[i], read, start+path.str.size()-kgraph,numbers)){
						return true;
					}
				}
				numbers.resize(size);
			}
		}
	}
	return false;
}


bool MappingSupervisor::alignOnPathFathers2(const path& path, const string& read, size_t position, vector<uNumber>& numbers ){
	//		cout<<"father : "<<path.str<<endl;
	//	if(read.size()<=position+offset){return true;}
	int start ((int)read.size()-(int)position-(int)path.str.size()*1);

	string region(read.substr(max(0,start),path.str.size()*2));
	//	if(region.size()<offset){return true;}

	//		unordered_set<minimizer> genomicKmers=allKmerSet(k2,homocompression(path.str));
	unordered_set<minimizer> genomicKmers=allKmerSet(k2,path.str);
	//	cout<<region<<" "<<path.str<<endl;
	//		if(jaccard(k2,homocompression(region),genomicKmers)>=minJacc){
	if(jaccard(k2,region,genomicKmers)>=minJacc){
		//		cout<<"good"<<endl;
		regionmapped++;
		if(read.size()<=position+path.str.size()-kgraph+offset){
			return true;
		}else{
			size_t size(numbers.size());
			numbers.insert(numbers.end(),path.numbers.begin(),path.numbers.end());
			auto list(listPathFathers(offset, path.str.substr(0,kgraph),0));
			for(size_t i(0);i<list.size();++i){
				if(alignOnPathFathers(list[i], read, position+path.str.size()-kgraph,numbers)){
					return true;
				}
			}
			numbers.resize(size);
		}
	}
	return false;
}

bool MappingSupervisor::alignOnPathFathers(const path& path, const string& read, size_t position, vector<uNumber>& numbers ){
	unordered_set<minimizer> genomicKmers=allKmerSet(k2,path.str);
	minimizer kmer;
	int start(0);
	bool found(false);
	for(uint i(0); i<path.str.size(); ++i){
		if(read.size()>k2+position+i){
			kmer=seq2int(read.substr(read.size()-k2-position-i,k2));
			if(genomicKmers.unordered_set::count(kmer)!=0){
				start=(int)read.size()-(int)k2-(int)position-i;
				found=true;
				break;
			}
		}
	}

	if(found){
		string region(read.substr(max(0,start),path.str.size()));
		if(jaccard(k2,region,genomicKmers)>=minJacc){
			//		cout<<"good"<<endl;
			regionmapped++;
			if(read.size()<=position+path.str.size()-kgraph+offset){
				return true;
			}else{
				size_t size(numbers.size());
				numbers.insert(numbers.end(),path.numbers.begin(),path.numbers.end());
				auto list(listPathFathers(offset, path.str.substr(0,kgraph),0));
				for(size_t i(0);i<list.size();++i){
					if(alignOnPathFathers(list[i], read, start+path.str.size()-kgraph,numbers)){
						return true;
					}
				}
				numbers.resize(size);
			}
		}
	}
	return false;
}




void MappingSupervisor::MapPart(size_t L, size_t R){
	//foreach unitig (sort of)
	for (size_t i(L); i<R; ++i){
		string unitig=unitigs[i];
		//		cout<<"unitig :"<<unitig<<endl;
		//		cin.get();
		if(unitig.size()>=minSizeUnitigs){
			if(bigUnitig++%100==0){
				cout<<bigUnitig++<<endl;
				cout<<"Unitigs mapped "<<unitigsPreMapped<<" Percent unitigs mapped : "<<(100*unitigsPreMapped)/(bigUnitig+1)<<endl;
				cout<<"Read pre-mapped : "<<readMapped<<" Percent read pre-mapped : "<<(100*readMapped)/(reads.size()+1)<<endl;
				cout<<"Read mapped on graph: "<<aligneOnPathSucess<<" Percent read mapped  on graph: "<<(100*aligneOnPathSucess)/(reads.size()+1)<<endl;
				cout<<island<<" islands..."<<endl;
				cout<<regionmapped*offset/(1000*1000)<<" Mnuc mapped"<<endl;
				cout<<leftmap<<" "<< rightmap<<endl;
				cout<<leftmapFail<<" "<<rightmapFail<<endl;
				cout<<endl;
			}
			vector<path> ListSons(listPathSons(offset, unitigs[i].substr(unitigs[i].size()-kgraph,kgraph),0));
			vector<path> ListFathers(listPathFathers(offset, unitigs[i].substr(0,kgraph),0));
			if( ListFathers.size()==0 and ListSons.size()==0){
				//								cin.get();
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
								//cout<<path.str<<endl;
								if(alignOnPathSons(path, read, 0,ref) or alignOnPathFathers(path, read, 0,ref)){
									mappedOnGraph=true;
								}

								if(mappedOnGraph){
									++aligneOnPathSucess;
									//									cout<<unitig<<endl;
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
						//						cout<<"try candidate correct"<<endl;
						int position,position1,position2;
						bool correct1(isCandidateCorrect(unitig,it->first,read2min,genomicKmers,position1));
						bool correct2(isCandidateCorrectReverse(unitig,it->first,read2min,genomicKmersRC,position2));
						string read;
						if(correct1){
							position=max(0,position1);
							read=(reads[it->first]);
							if(read.empty())
								continue;
						}else{
							if(correct2) {
								position=max(position2,0);
								read=reversecomplement(reads[it->first]);
								if(read.empty())
									continue;
							}else{
								continue;
							}
						}
						//						cout<<position<<endl;
						//						cout<<read.substr(0,position+kgraph)<<endl;
						//						cout<<read.substr(position+unitig.size()-kgraph)<<endl;
						//the read is mapped on the unitig (pre-mapped)
						if(!done){
							unitigsPreMapped++;
							done=true;
						}
						bool mappedLeft(false),mappedRight(false);

						vector<uNumber> ref,ref2;
						//For each possible path
						if(position+kgraph<offset){
							mappedLeft=true;
							//							leftmap++;
						}else{
							for(size_t j(0); j<ListFathers.size(); ++j){
								path pathFather(ListFathers[j]);
								//								cout<<pathFather.str<<endl;
								if(alignOnPathFathers(pathFather, read.substr(0,position+kgraph), 0,ref)){
									mappedLeft=true;
									leftmap++;
									break;
								}
							}
							leftmapFail++;
						}
						if(read.size()<position+unitig.size()-kgraph+offset){
							mappedRight=true;
							//							cout<<"rightmap"<<endl;
							//							rightmap++;
						}else{
							//							cout<<ListSons.size()<<endl;
							for(size_t j(0); j<ListSons.size(); ++j){
								path pathSon(ListSons[j]);
								if(alignOnPathSons(pathSon, read.substr(position+unitig.size()-kgraph), 0,ref2)){
									mappedRight=true;
									rightmap++;
									break;
								}
							}
							rightmapFail++;
						}

						if(mappedRight and mappedLeft){
							++aligneOnPathSucess;
							//							cout<<"success"<<endl;
							//							cin.get();
							//							cout<<"go2"<<endl;
							//							cout<<getPathBegin(ref)<<endl;
							//							cout<<getPathEnd(ref2)<<endl;
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



string MappingSupervisor::getPathEnd(vector<uNumber>& numbers){
	if(numbers.empty()){
		return "";
	}
	string path(unitigs[numbers[0]]);

	for(size_t i(1); i<numbers.size(); ++i){
		//		cout<<path<<endl;
		string unitig(unitigs[numbers[i]]);

		string inter=compactionEnd(path, unitig, kgraph);
		if(inter.empty()){
			path=compactionEnd(reversecomplement(path), unitig, kgraph);
		}else{
			path=inter;
		}
	}
	return path;
}

string MappingSupervisor::getPathBegin(vector<uNumber>& numbers){
	if(numbers.empty()){
		return "";
	}
	string path(unitigs[numbers[0]]);

	for(size_t i(1); i<numbers.size(); ++i){
		//		cout<<path<<endl;
		string unitig(unitigs[numbers[i]]);

		string inter=compactionBegin(path, unitig, kgraph);
		if(inter.empty()){
			path=compactionBegin(reversecomplement(path), unitig, kgraph);
		}else{
			path=inter;
		}
	}
	return path;
}





void MappingSupervisor::MapAll(){
	size_t nbThreads(4);
	vector<thread> threads;
	vector<size_t> limits = bounds(nbThreads, unitigs.size());

	for(size_t i(0); i<nbThreads; ++i) {
		threads.push_back(thread(&MappingSupervisor::MapPart, this, limits[i], limits[i+1]));
	}

	for(auto &t : threads){t.join();}
	cout<<"Unitigs mapped "<<unitigsPreMapped<<" Percent unitigs mapped : "<<(100*unitigsPreMapped)/(bigUnitig+1)<<endl;
	cout<<"Read pre-mapped : "<<readMapped<<" Percent read pre-mapped : "<<(100*readMapped)/(reads.size()+1)<<endl;
	cout<<"Read mapped on graph: "<<aligneOnPathSucess<<" Percent read mapped  on graph: "<<(100*aligneOnPathSucess)/(reads.size()+1)<<endl;
	cout<<island<<" islands..."<<endl;
}





