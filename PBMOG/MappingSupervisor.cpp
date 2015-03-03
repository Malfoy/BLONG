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
	if(depth>depthMax){return {};}
	vector<path> paths;
	vector<uNumber> Sons(G.getBegin(substr));
	for (size_t i(0); i<Sons.size(); ++i){
		string son(unitigs[Sons[i]]);
		string newson(compactionEnd(substr, son, kgraph));
		if(newson.empty()){
			cout<<"wtfnewson"<<endl;
			cin.get();
		}
		if (son.size()>=lengthRequired){
			paths.push_back(path{newson,{Sons[i]}});
		}else{
			vector<path> recurPaths(listPathSons(lengthRequired-son.size()+kgraph, newson.substr(newson.size()-kgraph), depth+1));
			for (size_t j(0); j<recurPaths.size(); ++j){
				string sonAndCo(compactionEnd(newson,recurPaths[j].str, kgraph));
				if(sonAndCo.size()>=lengthRequired){
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

	for (size_t i(0); i<Fathers.size(); ++i){
		string father(unitigs[Fathers[i]]);
		string newfather(compactionBegin(substr, father, kgraph));
		if(newfather.empty()){
			cout<<"wtfnewfather"<<endl;
			cin.get();
		}

		if (father.size()>=lengthRequired){
			paths.push_back(path{newfather,{Fathers[i]}});
		}else{
			vector<path> pathsRecur(listPathFathers(lengthRequired-father.size()+kgraph,newfather.substr(0,kgraph), depth+1));
			for (size_t j(0); j<pathsRecur.size(); ++j){
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


void MappingSupervisor::findCandidate(const string& unitig, unordered_set<minimizer>& minSet, unordered_map<rNumber,size_t>& Candidate, unordered_map<rNumber,unordered_set<minimizer>>& read2Min){
	if(unitigs.size()<H){
		minimizer kmerS=seq2intStranded(unitig.substr(0,k));
		minimizer kmerRC=seq2intStranded(reversecomplement(unitig.substr(0,k)));
		minimizer seq(min(kmerRC,kmerS));
		for(size_t i(0);;++i){
			if(minSet.unordered_set::count(seq)==0){
				minSet.insert(seq);
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
			if(i+k<unitig.size()){
				updateMinimizer(kmerS, unitig[i+k], k);
				updateMinimizerRC(kmerRC, unitig[i+k], k);
				seq=(min(kmerRC,kmerS));
			}else{
				return;
			}
		}
	}else{
		vector<minimizer> sketch(minHashpart(H,k,unitig,part));
		//For each minimizer of the unitig
		for(size_t j(0);j<sketch.size();++j){
			minimizer seq=sketch[j];
			if(minSet.unordered_set::count(seq)==0){
				minSet.insert(seq);
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



bool MappingSupervisor::isCandidateCorrect(const string& unitig, rNumber readNumber, unordered_map<rNumber,unordered_set<minimizer>>& read2min, unordered_set<minimizer>& genomicKmers,int& position,bool strand){
	unordered_set<minimizer> setMin(read2min[readNumber]);
	string read;
	if(strand){
		read=reads[readNumber];
	}else{
		read=reversecomplement(reads[readNumber]);
	}
	if(read.empty()){return false;}
	minimizer minS(seq2intStranded(read.substr(0,k)));
	minimizer minRC(seq2intStranded(reversecomplement(read.substr(0,k))));
	minimizer mini(min(minS,minRC));
	for(size_t j(0);;++j){
		if(setMin.unordered_set::count(mini)!=0){
			int posInSeq(positionInSeq(unitig, mini, k));
			if(posInSeq>=0){
				int pos((int)j-posInSeq);
				string region(read.substr(max(pos,0),1*unitig.size()));
				if(genomicKmers.empty()){
					genomicKmers=allKmerSetStranded(k2,unitig);
				}
				if(jaccardStranded(k2,region,genomicKmers)>=minJacc){
					readMapped++;
					position=(pos);
					return true;

				}
			}
		}
		if(j+k<read.size()){
			updateMinimizer(minS, read[j+k], k);
			updateMinimizerRC(minRC, read[j+k], k);
			mini=min(minS,minRC);
		}else{
			return false;
		}
	}
	return false;
}






bool MappingSupervisor::alignOnPathSons(const path& path, const string& read, size_t position,vector<uNumber>& numbers){
	unordered_set<minimizer> genomicKmers=allKmerSetStranded(k2,path.str);
//	cout<<"go"<<endl;
	minimizer kmer(seq2intStranded(read.substr(position,k2)));
//	cout<<"kmer done"<<endl;
	int start(0);
	bool found(false);

	for(uint i(0); i<path.str.size(); ++i){
		if(genomicKmers.unordered_set::count(kmer)!=0){
			start=max(((int)i+(int)position-(int)positionInSeqStranded(path.str, kmer, k2)),0);
			found=true;
			break;
		}
		if(read.size()>position+i+k2){
			updateMinimizer(kmer, read[position+i+k2], k2);
		}else{
			return false;
		}
	}

	if(found){
//		cout<<"found"<<endl;
		string region(read.substr(start,path.str.size()));
//		cout<<region<<endl;
//		cout<<"start ok"<<endl;
		if(jaccardStranded(k2,region,genomicKmers)>=minJacc){
			regionmapped++;
			size_t size(numbers.size());
			numbers.insert(numbers.end(),path.numbers.begin(),path.numbers.end());
			if(read.size()<start+path.str.size()-kgraph+offset){
				return true;
			}else{
//				cout<<read.size()<<endl;
//				cout<<start+path.str.size()-kgraph<<endl;
//				cout<<path.str<<endl;
				auto list(listPathSons(offset, path.str.substr(path.str.size()-kgraph,kgraph),0));
//				cout<<"recur"<<endl;
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


bool MappingSupervisor::alignOnPathFathers(const path& path, const string& read, size_t position, vector<uNumber>& numbers ){
	unordered_set<minimizer> genomicKmers=allKmerSetStranded(k2,path.str);
	minimizer kmer(seq2intStranded(read.substr(read.size()-k2-position,k2)));
	int start(0);
	bool found(false);
	for(uint i(0); i<path.str.size(); ++i){
		if(genomicKmers.unordered_set::count(kmer)!=0){
			start=max(((int)read.size()-(int)position-(int)i-(int)positionInSeqStranded(path.str, kmer, k2)),0);
			found=true;
			break;
		}
		if(read.size()>k2+position+i){
			updateMinimizer(kmer, read[read.size()-k2-position-i], k2);
		}else{
			return false;
		}
	}

	if(found){
		string region(read.substr(start,path.str.size()));
		if(jaccardStranded(k2,region,genomicKmers)>=minJacc){
			regionmapped++;
			size_t size(numbers.size());
			numbers.insert(numbers.end(),path.numbers.begin(),path.numbers.end());
			if(read.size()<start+path.str.size()-kgraph+offset){
				return true;
			}else{
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
			if( ListFathers.empty() and ListSons.empty()){
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
					if(it->second>=multi){
						candidate++;
						int position,position1,position2;
						bool correct1(isCandidateCorrect(unitig,it->first,read2min,genomicKmers,position1,true));
						string read;
						if(correct1){
							position=max(0,position1);
							read=(reads[it->first]);
							if(read.empty()){
								continue;
							}
						}else{
							bool correct2(isCandidateCorrect(unitig,it->first,read2min,genomicKmersRC,position2,false));
							if(correct2) {
								position=max(position2,0);
								read=reversecomplement(reads[it->first]);
								if(read.empty()){
									continue;
								}
							}else{
								continue;
							}
						}

						//the read is mapped on the unitig (pre-mapped)
						if(!done){
							unitigsPreMapped++;
							done=true;
						}

						//						continue;

						bool mappedLeft(false),mappedRight(false);

						vector<uNumber> ref,ref2;
						//For each possible path
						if(position+kgraph<offset){
							mappedLeft=true;
						}else{
							for(size_t j(0); j<ListFathers.size(); ++j){
								path pathFather(ListFathers[j]);
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
						}else{
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
		string unitig(unitigs[numbers[i]]);

		string inter=compactionEnd(path, unitig, kgraph);
		if(inter.empty()){
			path=compactionEnd(path, reversecomplement(unitig), kgraph);
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
		string unitig(unitigs[numbers[i]]);
		string inter=compactionBegin(path, unitig, kgraph);
		if(inter.empty()){
			path=compactionBegin(path, reversecomplement(unitig), kgraph);
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
	cout<<candidate<<" candidate..."<<endl;
}





