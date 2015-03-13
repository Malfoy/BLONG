//
//  MappingSupervisor.cpp
//  PBMOG

//
//  Created by malfoy on 09/01/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#include "MappingSupervisor.h"
#include "Utils.h"
#include "nw.h"






vector<path> MappingSupervisor::listPathSons(size_t lengthRequired, const string& substr, int depth){
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



vector<path> MappingSupervisor::listPathFathers(size_t lengthRequired, const string& substr, int depth){
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


void MappingSupervisor::findCandidate(const string& unitig, unordered_set<minimizer>& minSet, unordered_map<rNumber,size_t>& Candidate, vector<unordered_set<minimizer>>& read2Min){
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



bool MappingSupervisor::isCandidateCorrect2(const string& unitig, const string& read, unordered_set<minimizer>& genomicKmers,int& position, unordered_set<minimizer>& setMin){
	if(read.empty()){return false;}
	minimizer minS(seq2intStranded(read.substr(0,k)));
	minimizer minRC(seq2intStranded(reversecomplement(read.substr(0,k))));
	minimizer mini(min(minS,minRC));
	for(int i(0);;++i){
		if(setMin.unordered_set::count(mini)!=0){
			int posInSeq(positionInSeq(unitig, mini, k));
			if(posInSeq>=0){

				int pos(max((i-posInSeq),0));
				string region(read.substr(pos,unitig.size()));
				if(jaccardStranded(k2,region,genomicKmers)>=minJacc){
					readMapped++;
					position=(pos);
					return true;
				}
			}
		}
		if(i+k<read.size()){
			updateMinimizer(minS, read[i+k], k);
			updateMinimizerRC(minRC, read[i+k], k);
			mini=min(minS,minRC);
		}else{
			return false;
		}
	}
	return false;
}


bool MappingSupervisor::isCandidateCorrect(const string& unitig, const string& read, unordered_set<minimizer>& genomicKmers,int& position, unordered_set<minimizer>& setMin){
	if(read.empty()){return false;}
	minimizer minS(seq2intStranded(read.substr(0,k2)));
	//	minimizer minRC(seq2intStranded(reversecomplement(read.substr(0,k))));
	//	minimizer mini(min(minS,minRC));
	for(int i(0);;++i){
		if(genomicKmers.unordered_set::count(minS)!=0){
			int posInSeq(positionInSeqStranded(unitig, minS, k2));
			if(posInSeq>=0){

				int pos(i-posInSeq);
				string region;
				if(pos>=0){
					region=(read.substr(pos,unitig.size()));
				}else{
					region=(read.substr(0,unitig.size()+pos));
				}
				//				cout<<"jacc "<<jaccardStranded(k2,region,genomicKmers)<<endl;
				//				cout<<region<<endl<<unitig<<endl;
				//				cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
				//				cout<<pos<<" "<<posInSeq<<" "<<i<<endl;
				//				cin.get();
				if(jaccardStranded(k2,region,genomicKmers)>=minJacc){

					//					cout<<"candidate correct"<<endl;
					//					cout<<region<<" "<<unitig<<endl;
					//					cin.get();
					readMapped++;
					position=(pos);
					return true;
				}
			}
		}
		if(i+k2<read.size()){
			updateMinimizer(minS, read[i+k2], k2);
			//			updateMinimizerRC(minRC, read[i+k], k);
			//			mini=min(minS,minRC);
		}else{
			return false;
		}
	}
	return false;
}





bool MappingSupervisor::alignOnPathSons(const path& path, const string& read, size_t position,vector<uNumber>& numbers){
	unordered_set<minimizer> genomicKmers=allKmerSetStranded(k2,path.str);
	minimizer kmer(seq2intStranded(read.substr(position,k2)));
	int start(0);
	bool found(false);

	for(uint i(0); i<path.str.size();++i){
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
		string region(read.substr(start,path.str.size()));
		//		cout<<"region : "<<region<<" path : "<<path.str<<endl;
		if(jaccardStranded(k2,region,genomicKmers)>=minJacc){
			size_t size(numbers.size());
			for(int i((int)path.numbers.size()-1);i>=0;--i){
				numbers.push_back(path.numbers[i]);
			}
			//			numbers.insert(numbers.end(),path.numbers.begin(),path.numbers.end());
			if(read.size()<start+path.str.size()+offset){
				return true;
			}else{
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


bool MappingSupervisor::alignOnPathFathers(const path& path, const string& read, size_t position, vector<uNumber>& numbers ){
	//	cout<<read<<endl;
	unordered_set<minimizer> genomicKmers=allKmerSetStranded(k2,path.str);
	minimizer kmer(seq2intStranded(read.substr(position-k2,k2)));
	int start(0);
	bool found(false);
	for(uint i(0); i<path.str.size();){
		//		printMinimizer(kmer, k2);
		if(genomicKmers.unordered_set::count(kmer)!=0){
			start=max(((int)position-(int)i-positionInSeqStrandedEnd(path.str, kmer, k2)),0);
			found=true;
			break;
		}
		++i;
		if(position>=k2+i){
			updateMinimizerEnd(kmer, read[position-k2-i], k2);
		}else{
			return false;
		}
	}

	if(found){
		string region(read.substr(max(start-(int)k2,0),path.str.size()));
		if(jaccardStranded(k2,region,genomicKmers)>=minJacc){
			//			regionmapped+=region.size();
			size_t size(numbers.size());
			numbers.insert(numbers.end(),path.numbers.begin(),path.numbers.end());
			//			path.numbers.insert(path.numbers.end(),numbers.begin(),numbers.end());
			if(start<(int)offset){
				return true;
			}else{
				auto list(listPathFathers(offset, path.str.substr(0,kgraph),0));
				for(size_t i(0);i<list.size();++i){
					if(alignOnPathFathers(list[i], read, start,numbers)){
						return true;
					}
				}
				numbers.resize(size);
			}
		}
	}
	return false;
}

void MappingSupervisor::MapFromUnitigs(const string& unitig){
	unordered_set<minimizer> min;
	unordered_map<rNumber,size_t> Candidate;Candidate.set_empty_key(-1);
	vector<unordered_set<minimizer>> read2min(reads.size());
	findCandidate(unitig,min,Candidate,read2min);
	bool done(false);
	unordered_set<minimizer> genomicKmers(allKmerSetStranded(k2,unitig));
	vector<path> ListSons(listPathSons(offset, unitig.substr(unitig.size()-kgraph,kgraph),0));
	//	vector<path> ListFathers(listPathFathers(offset, unitig.substr(0,kgraph),0));
	vector<path> ListFathers(listPathSons(offset, reversecomplement(unitig.substr(0,kgraph)),0));
	if( ListFathers.empty() and ListSons.empty()){
		island++;
		return;
	}
	if(Candidate.empty()){
		return;
	}
	//foreach reads that could map on the unitig (prepremapped)
	for(auto it=Candidate.begin(); it!=Candidate.end(); ++it){

		if(it->second>=multi){
			//						readused.insert(it->first);
			string read=reads[it->first];
			if(read.empty()){
				continue;
			}
			bool preMapped(false);
			bool stranded;
			int position,position1(-1),position2(-1);

			bool correct1(isCandidateCorrect(unitig,read,genomicKmers,position1,read2min[it->first]));
			if(correct1){
				position=max(0,position1);
				preMapped=true;
				stranded=true;
			}else{
				read=reversecomplement(read);
				bool correct2(isCandidateCorrect(unitig,read,genomicKmers,position2,read2min[it->first]));
				if(correct2){
					position=max(position2,0);
					preMapped=true;
					stranded=false;
				}else{
					candidate++;
					continue;
				}
			}

			if(preMapped){
				if(!done){
					done=true;
					//regionmapped+=unitig.size();
					unitigsPreMapped++;
				}
				bool mappedLeft(false),mappedRight(false);

				vector<uNumber> numberBegin,numberEnd;
				string beg,end;
				if(position+kgraph<offset){
					mappedLeft=true;
					//leftmap++;
				}else{
					beg=(read.substr(0,position+kgraph));
				}
				if(read.size()<position+unitig.size()+offset-kgraph){
					mappedRight=true;
					//	rightmap++;
				}else{
					end=(read.substr(position+unitig.size()-kgraph));
				}
				if(!mappedLeft){
					//					cout<<read<<endl;
					//					cout<<"beg : "<<reversecomplement(beg)<<endl;
					//					cout<<ListFathers.size()<<endl;
					for(size_t j(0); j<ListFathers.size(); ++j){
						path pathFather(ListFathers[j]);
						//						pathFather.str=reversecomplement(pathFather.str);
						//						cout<<"pather : "<<pathFather.str<<endl;
						//						cout<<reversecomplement(beg)<<endl;
						//						cout<<"pather : "<<reversecomplement(pathFather.str)<<endl;
						if(alignOnPathSons(pathFather, reversecomplement(beg), 0,numberBegin)){
							//							cout<<"sucess"<<endl;
							//							cout<<"sucess"<<endl;cin.get();
							//regionmapped+=beg.size();
							mappedLeft=true;
							if(mapPartAllowed){
								mutexEraseReads.lock();
								reads[it->first]=end;
								mutexEraseReads.unlock();
							}
							leftmap++;
							break;
						}

					}
					leftmapFail++;
					//					cout<<"FAIl"<<endl;cin.get();
				}
				if(!mappedRight){
					for(size_t j(0); j<ListSons.size(); ++j){
						path pathSon(ListSons[j]);
						//						cout<<"pathson : "<<pathSon.str<<endl;
						if(alignOnPathSons(pathSon,end , 0,numberEnd)){
							//							cout<<"sucess"<<endl;cin.get();
							//regionmapped+=end.size();
							mappedRight=true;
							if(mapPartAllowed){
								mutexEraseReads.lock();
								reads[it->first]=beg;
								mutexEraseReads.unlock();
							}
							rightmap++;
							break;
						}

					}
					rightmapFail++;
					//					cout<<"FAIl"<<endl;cin.get();
				}

				if(mappedRight and mappedLeft){
//					cout<<"go "<<unitig<<endl;
//					cout<<numberBegin.size()<<" "<<numberEnd.size()<<endl;
					string pathBegin(getPathEnd(numberBegin));
//					cout<<pathBegin<<endl;
					pathBegin=pathBegin.substr(0,beg.size());
//					cout<<pathBegin<<endl;
					if(pathBegin.empty() and numberBegin.size()!=0){
						cout<<"fail to recompose path begin"<<endl;
					}

					string pathEnd(getPathEnd(numberEnd));
//					cout<<pathEnd<<endl;
					pathEnd=pathEnd.substr(0,end.size());
//					cout<<pathEnd<<endl;
					if(pathEnd.empty() and numberEnd.size()!=0){
						cout<<"fail to recompose path end"<<endl;
					}

					string path;
					if(!numberBegin.empty()){
						path=(compactionBegin(unitig,pathBegin,kgraph));
						if(path.empty()){
							cout<<"fail compactbegin..."<<endl;
							path=unitig;
						}
					}else{
//						path=unitig.substr(max(0,(int)unitig.size()-((int)read.size()-(int)end.size())));
						if(stranded){
							if(position1>0){
								path=unitig;
							}else{
								path=unitig.substr(-position1);
							}
						}else{
							if(position2>0){
								path=unitig;
							}else{
								path=unitig.substr(-position2);
							}
						}
					}

					if(!numberEnd.empty()){
						string final(compactionEnd(path, pathEnd, kgraph));
						if(!final.empty()){
							path=final;
						}else{
							cout<<"fail compactend..."<<endl;
						}
					}else{
						path=path.substr(0,read.size());
					}


					if(path.empty()){
						cout<<"Cant recompose path...."<<endl;
					}


					string seq1,seq2;
//					cout<<"nw"<<endl;
//					cout<<"!"<<read<<"!"<<" !"<<path<<"!"<<endl;
					nw(path, read, seq1, seq2, false);
					if(scoreFromAlignment(seq1)<=errorRate){
//						cout<<read<<" ENDREAD"<<endl;
//						cout<<path<<" ENDPATH"<<endl;
//						cout<<seq1<<" ENDALIGN"<<endl;
//						cout<<seq2<<" ENDALIGN"<<endl;
//						cin.get();
						mutexEraseReads.lock();
						if(!reads[it->first].empty()){
							outFile<<"read : "<<read<<endl;
							outFile<<"path : "<<path<<endl;
//							reads[it->first].clear();
							++aligneOnPathSucess;
							regionmapped+=read.size();
						}
						mutexEraseReads.unlock();
					}else{
						fail++;
						cout<<read<<" "<<path<<" "<<seq1<<endl;
//						cout<<scoreFromAlignment(seq1)<<endl;
//						cin.get();
					}

					continue;
				}
			}
		}
	}
}




void MappingSupervisor::MapPart(){
	//foreach unitig (sort of)
	while(indice<unitigs.size()){
		mutexReadReads.lock();
		const string unitig=unitigs[indice];
		indice++;
		mutexReadReads.unlock();
		//		cout<<"unitig :"<<unitig<<endl;
		if(unitig.size()>=minSizeUnitigs){
			bigUnitig++;
			if(bigUnitig%1==0){
				cout<<bigUnitig<<endl;
				cout<<"Unitigs mapped "<<unitigsPreMapped<<" Percent unitigs mapped : "<<(100*unitigsPreMapped)/(bigUnitig)<<endl;
				cout<<"Read pre-mapped : "<<readMapped<<" Percent read pre-mapped : "<<(100*readMapped)/(reads.size())<<endl;
				cout<<"Read used: "<<readMapped<<" Percent read used : "<<(100*readused.size())/(reads.size())<<endl;
				cout<<"Read mapped on graph: "<<aligneOnPathSucess<<" Percent read mapped  on graph: "<<(100*aligneOnPathSucess)/(reads.size())<<endl;
				cout<<island<<" islands..."<<endl;
				cout<<regionmapped/(1000*1000)<<" Mnuc mapped or "<<regionmapped/(1000*1)<<" Knuc "<<endl;
				cout<<leftmap<<" "<< rightmap<<endl;
				cout<<leftmapFail<<" "<<rightmapFail<<endl;
				if(fail+aligneOnPathSucess!=0){
					cout<<(100*aligneOnPathSucess)/(fail+aligneOnPathSucess);
				}
				cout<<endl;
			}

			if(unitig.size()<offset){

				vector<path> ListFathers(listPathFathers(offset, unitig.substr(0,kgraph),0));
				for(size_t j(0); j<ListFathers.size(); ++j){
					path pathFather(ListFathers[j]);
					string newUnitig(compactionBegin(unitig, pathFather.str, kgraph));
					//check if not wtf
					MapFromUnitigs(newUnitig);
				}
				vector<path> ListSons(listPathSons(offset, unitig.substr(unitig.size()-kgraph,kgraph),0));
				for(size_t j(0); j<ListSons.size(); ++j){
					path pathSon(ListSons[j]);
					string newUnitig(compactionEnd(unitig, pathSon.str, kgraph));
					//check if not wtf
					MapFromUnitigs(newUnitig);
				}

			}else{
				MapFromUnitigs( unitig);
			}
		}
	}
}



string MappingSupervisor::getPathEnd(vector<uNumber>& numbers){
	if(numbers.empty()){
		return "";
	}
	//	for(size_t i(0); i<numbers.size(); ++i){
	//		string unitig(unitigs[numbers[i]]);
	//		cout<<unitig<<endl;
	//	}
	//	cout<<"end"<<endl;

	string path(unitigs[numbers[0]]);

	for(size_t i(1); i<numbers.size(); ++i){
		string unitig(unitigs[numbers[i]]);

		string inter=compactionEnd(path, unitig, kgraph);
		if(inter.empty()){
			path=compactionEnd(reversecomplement(path), unitig, kgraph);
			//			if(path.empty()){
			//				cout<<"!"<<path<<"!"<<unitig<<endl;
			//				cin.get();
			//			}
		}else{
			path=inter;
		}
		if(path.empty()){
		}
	}
	return path;
}


string MappingSupervisor::getPathBegin(vector<uNumber>& numbers){
	if(numbers.empty()){
		return "";
	}
	//	reverse(numbers.begin(), numbers.end());
	string path(unitigs[numbers[0]]);

	for(size_t i(1); i<numbers.size(); ++i){
		string unitig(unitigs[numbers[i]]);
		//		cout<<unitig<<" : unitig "<<endl;
		string inter=compactionBegin(path, unitig, kgraph);
		if(inter.empty()){
			inter=compactionBegin(reversecomplement(path), unitig, kgraph);
			if(inter.empty()){
				cout<<"wtf : "<<path<<" "<<unitig<<endl;
				cin.get();
			}else{
				path=inter;
			}
		}else{
			path=inter;
		}
		//		cout<<inter<<" <-this was inter lol"<<endl;
	}
	return path;
}





void MappingSupervisor::MapAll(){
	vector<thread> threads;
	vector<size_t> limits = bounds(nbThreads, unitigs.size());

	for(size_t i(0); i<nbThreads; ++i) {
		threads.push_back(thread(&MappingSupervisor::MapPart, this));
	}

	for(auto &t : threads){t.join();}
	
	//	for(size_t i(0);i<reads.size();++i){
	//		if(!reads[i].empty()){
	//			cout<<reads[i]<<endl;
	//			cin.get();
	//		}
	//	}
	cout<<bigUnitig<<endl;
	cout<<"Unitigs mapped "<<unitigsPreMapped<<" Percent unitigs mapped : "<<(100*unitigsPreMapped)/(bigUnitig)<<endl;
	cout<<"Read pre-mapped : "<<readMapped<<" Percent read pre-mapped : "<<(100*readMapped)/(reads.size())<<endl;
	cout<<"Read used: "<<readMapped<<" Percent read used : "<<(100*readused.size())/(reads.size())<<endl;
	cout<<"Read mapped on graph: "<<aligneOnPathSucess<<" Percent read mapped  on graph: "<<(100*aligneOnPathSucess)/(reads.size())<<endl;
	cout<<island<<" islands..."<<endl;
	cout<<regionmapped/(1000*1000)<<" Mnuc mapped or "<<regionmapped/(1000*1)<<" Knuc "<<endl;
	cout<<leftmap<<" "<< rightmap<<endl;
	cout<<leftmapFail<<" "<<rightmapFail<<endl;
	cout<<endl;
}





