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
#include <math.h>






vector<path> MappingSupervisor::listPathSons(size_t lengthRequired, const string& substr, int depth){
	if(depth>depthMax){deepper++;return {};}
	vector<path> paths;
	vector<uNumber> Sons(G.getBegin(substr));
	for (size_t i(0); i<Sons.size(); ++i){
		string son(unitigs[Sons[i]]);
		//		cout<<son<<endl;
		string newson(compactionEnd(substr, son, kgraph));

		if(newson.empty()){cout<<"wtfnewson"<<endl;cin.get();}//DEBUG

		if (son.size()>=lengthRequired){
			paths.push_back(path{newson,{Sons[i]}});
		}else{
			vector<path> recurPaths(listPathSons(lengthRequired-son.size()+kgraph, newson.substr(newson.size()-kgraph), depth+1));
			for (size_t j(0); j<recurPaths.size(); ++j){
				string sonAndCo(compactionEnd(newson,recurPaths[j].str, kgraph));
				if(sonAndCo.size()>=lengthRequired){ //DEBUG
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


string MappingSupervisor::getRead(rNumber n){
	string res,more("");
	rPosition pos(number2position[n]);
	mutexReadReads.lock();
	reads.seekg(pos, ios::beg);
	getline(reads,res);
	while(reads.peek()!='>' and !reads.eof()){
		getline(reads,more);
		res+=more;
	}
	if(reads.eof()){reads.clear();}
	mutexReadReads.unlock();
	return res;
}


void MappingSupervisor::findCandidate(const string& unitig, unordered_map<rNumber,uint32_t>& Candidate, unordered_map<rNumber,unordered_set<minimizer>>& read2Min){
	if(unitigs.size()<H){
		unordered_set<minimizer> minSet;
		vector<rNumber> readsNumbers;
		minimizer kmerS=seq2intStranded(unitig.substr(0,k));
		minimizer kmerRC=seq2intStranded(reversecomplement(unitig.substr(0,k)));
		minimizer seq(min(kmerRC,kmerS));
		for(size_t i(0);;++i){
			if(minSet.unordered_set::count(seq)==0){
				minSet.insert(seq);
				if(min2Reads.unordered_map::count(seq)!=0){
					readsNumbers=min2Reads[seq];
					for(size_t j(0);j<readsNumbers.size();++j){
						if(number2position[readsNumbers[j]]!=0){
							++Candidate[readsNumbers[j]];
							read2Min[readsNumbers[j]].insert(seq);
						}else{
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
		vector<rNumber> readsNumbers;
		removeDuplicate(sketch);
		//For each minimizer of the unitig
		for(size_t i(0);i<sketch.size();++i){
			minimizer seq(sketch[i]);
			if(min2Reads.unordered_map::count(seq)!=0){
				readsNumbers=min2Reads[seq];
				for(size_t j(0);j<readsNumbers.size();++j){
					if(number2position[readsNumbers[j]]!=0){
						++Candidate[readsNumbers[j]];
						read2Min[readsNumbers[j]].insert(seq);
					}
				}
			}
		}
	}
}



bool MappingSupervisor::isCandidateCorrect(const string& unitig, const string& read, unordered_set<minimizer>& genomicKmers,int& position, unordered_set<minimizer>& setMin){
	if(read.empty()){return false;}
	minimizer minS(seq2intStranded(read.substr(0,k)));
	minimizer minRC(seq2intStranded(reversecomplement(read.substr(0,k))));
	minimizer mini(min(minS,minRC));
	string region;
	for(int i(0);;++i){
		if(setMin.unordered_set::count(mini)!=0){
			int posInSeq(positionInSeq(unitig, mini, k));
			if(posInSeq>=0){
				int pos(max((i-posInSeq),0));
				region=(read.substr(pos,unitig.size()));
				//				cout<<region.size()<<" "<<log2(region.size())/2<<" "<<100*(pow(1-0.1,log2(region.size())/2))<<endl;
				//				size_t kAlt((size_t)log2(region.size())/2+1);
				//				double minJaccAlt(100*(pow(1-0.1,kAlt)));
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


bool MappingSupervisor::isCandidateCorrect2(const string& unitig, const string& read, unordered_set<minimizer>& genomicKmers,int& position, unordered_set<minimizer>& setMin){
	if(read.empty()){return false;}
	minimizer minS(seq2intStranded(read.substr(0,k2)));
	for(int i(0);;++i){
		if(genomicKmers.unordered_set::count(minS)!=0){
			int posInSeq(positionInSeqStranded(unitig, minS, k2));
			if(posInSeq>=0){
				int pos(i-posInSeq);
				string region=(pos>=0 ? read.substr(pos,unitig.size()): read.substr(0,unitig.size()+pos));
				if(jaccardStranded(k2,region,genomicKmers)>=minJacc){
					readMapped++;
					position=(pos);
					return true;
				}
			}
		}
		if(i+k2<read.size()){
			updateMinimizer(minS, read[i+k2], k2);
		}else{
			return false;
		}
	}
	return false;
}


//Can do way better much redundancy
bool MappingSupervisor::isCandidateCorrectMap(const string& unitig, const string& read, const unordered_multimap<string,string>& genomicKmers,int& position, unordered_set<minimizer>& setMin, int&positionRead){
	if(read.empty()){return false;}
	minimizer minS(seq2intStranded(read.substr(0,k)));
	minimizer minRC(seq2intStranded(reversecomplement(read.substr(0,k))));
	minimizer mini(min(minS,minRC));
	string region;

	for(int i(0);;++i){
		if(setMin.unordered_set::count(mini)!=0){
			int posInSeq(positionInSeq(unitig, mini, k));
			if(posInSeq>=0){
				int pos(i-posInSeq);
				region=(pos>=0 ? read.substr(pos,unitig.size()): read.substr(0,unitig.size()+pos));
				if(jaccardStrandedErrors(k2,region,genomicKmers,nuc)>=minJacc){
					readMapped++;
					position=(pos);
					positionRead=pos;
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
		if(jaccardStranded(k2,region,genomicKmers)>=minJacc){
			size_t size(numbers.size());
			for(int i((int)path.numbers.size()-1);i>=0;--i){
				numbers.push_back(path.numbers[i]);
			}
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


bool MappingSupervisor::alignOnPathsSons(const vector<path>& Paths, const string& read, size_t position,vector<uNumber>& numbers){
	double maxScore(0);
	size_t maxIndice(0);
	int maxStart(0);
	for(size_t ii(0); ii<Paths.size();++ii){
		path path(Paths[ii]);
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
				break;
			}
		}
		if(found){
			string region(read.substr(start,path.str.size()));
			double score(jaccardStranded(k2,region,genomicKmers));
			if(score>maxScore){
				maxScore=score;
				maxIndice=ii;
				maxStart=start;
			}
		}
	}
	if(maxScore>minJacc){
		path path(Paths[maxIndice]);
		for(int i((int)path.numbers.size()-1);i>=0;--i){
			numbers.push_back(path.numbers[i]);
		}
		return (read.size()<maxStart+path.str.size()+offset ? true : alignOnPathsSons(listPathSons(offset, path.str.substr(path.str.size()-kgraph,kgraph),0),read,maxStart+path.str.size()-kgraph,numbers));
	}

	return false;
}


bool MappingSupervisor::alignOnPathSonsErrors(const path& path, const string& read, size_t position,vector<uNumber>& numbers){
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
	auto genomicKmersErrors(allKmerMapStranded(k2,path.str,nuc));
	if(found){
		string region(read.substr(start,path.str.size()));
		if(jaccardStrandedErrors(k2,region,genomicKmersErrors,nuc)>=minJacc){
			size_t size(numbers.size());
			for(int i((int)path.numbers.size()-1);i>=0;--i){
				numbers.push_back(path.numbers[i]);
			}
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
			return false;
		}
	}
	return false;
}


bool MappingSupervisor::alignOnPathsSonsErrors(const vector<path>& Paths, const string& read, size_t position,vector<uNumber>& numbers){
	double maxScore(0);
	size_t maxIndice(0);
	int maxStart(0);
	for(size_t ii(0); ii<Paths.size();++ii){
		path path(Paths[ii]);
		unordered_set<minimizer> genomicKmers=allKmerSetStranded(k2,path.str);
		minimizer kmer(seq2intStranded(read.substr(position,k2)));
		int start(0);
		bool found(false);

		for(uint i(0); i<path.str.size();++i){
			if(genomicKmers.unordered_set::count(kmer)!=0){
				start=max(((int)i+(int)position-(int)positionInSeqStranded(path.str, kmer, k2)),0);
				if (start+(int)(path.str.size())-kgraph>position) {
					found=true;
					break;
				}
			}
			if(read.size()>position+i+k2){
				updateMinimizer(kmer, read[position+i+k2], k2);
			}else{
				break;
			}
		}
		if(found){
			auto genomicKmersErrors(allKmerMapStranded(k2,path.str,nuc));
			string region(read.substr(start,path.str.size()));
			double score(jaccardStrandedErrors(k2,region,genomicKmersErrors,nuc));
			if(score>maxScore){
				maxScore=score;
				maxIndice=ii;
				maxStart=start;
			}
		}
	}

	if(maxScore>minJacc){
		path path(Paths[maxIndice]);
		for(int i((int)path.numbers.size()-1);i>=0;--i){
			numbers.push_back(path.numbers[i]);
		}
		if (read.size()<maxStart+(int)(path.str.size())+offset) {
			return true;
		}else{
			if (alignOnPathsSonsErrors(listPathSons(offset, path.str.substr(path.str.size()-kgraph,kgraph),0),read,maxStart+(int)(path.str.size())-kgraph,numbers)) {
				return true;
			}else{
				return(maxStart+(path.str.size())>read.size());
			}
		}
	}
	return false;
}


bool MappingSupervisor::alignOnPathsSonsErrorsAll(const vector<path>& Paths, const string& read, size_t position,vector<uNumber>& numbers){
	for(size_t ii(0); ii<Paths.size();++ii){
		path path(Paths[ii]);
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
				break;
			}
		}

		if(found){
			auto genomicKmersErrors(allKmerMapStranded(k2,path.str,nuc));
			string region(read.substr(start,path.str.size()));
			double score(jaccardStrandedErrors(k2,region,genomicKmersErrors,nuc));
			if(score>minJacc){
				size_t size(numbers.size());
				for(int i((int)path.numbers.size()-1);i>=0;--i){
					numbers.push_back(path.numbers[i]);
				}
				if(read.size()<start+path.str.size()+offset){
					return true;
				}else{
					if(alignOnPathsSonsErrorsAll(listPathSons(offset, path.str.substr(path.str.size()-kgraph,kgraph),0),read,start+path.str.size()-kgraph,numbers)){
						return true;
					}
				}
				numbers.resize(size);
			}
		}
	}
	return false;
}


bool MappingSupervisor::preMapUnitig(const string& unitig, string& read,const unordered_multimap<string,string>& genomicKmers,int& position,unordered_set<minimizer> setMin, int& positionUnitig, bool& stranded, int& position1, int& position2){
	bool correct1(isCandidateCorrectMap(unitig,read,genomicKmers,position1,setMin,positionUnitig));
	if(correct1){
		position=max(0,position1);
		stranded=true;
		return true;
	}else{
		read=reversecomplement(read);
		bool correct2(isCandidateCorrectMap(unitig,read,genomicKmers,position2,setMin,positionUnitig));
		if(correct2){
			position=max(position2,0);
			stranded=false;
			return true;
		}else{
			candidate++;
		}
	}
	return false;
}


bool MappingSupervisor::mapOnGraph(vector<uNumber>& numberBegin, vector<uNumber>& numberEnd,const string& unitig, const string& read, int position, vector<path> ListFathers, vector<path> ListSons,rNumber rNum, string& beg, string& end){
	bool mappedLeft(false),mappedRight(false);
	if(position+kgraph<offset){
		mappedLeft=true;
	}else{
		beg=(read.substr(0,position+kgraph));
	}
	if(read.size()<position+unitig.size()+offset-kgraph){
		mappedRight=true;
	}else{
		end=(read.substr(position+unitig.size()-kgraph));
	}
	if(!mappedLeft){
		if(alignOnPathsSonsErrors(ListFathers, reversecomplement(beg), 0,numberBegin)){
			mappedLeft=true;
			leftmap++;
		}else{
			leftmapFail++;
		}
	}
	if(!mappedRight){
		if(alignOnPathsSonsErrors(ListSons,end , 0,numberEnd)){
			mappedRight=true;
			rightmap++;
		}else{
			rightmapFail++;

		}
	}
	//	return true;//TOREMOVE
	return mappedLeft and mappedRight;
}



string MappingSupervisor::recoverPath(vector<uNumber>& numberBegin, vector<uNumber>& numberEnd,size_t begsize,size_t endsize,size_t readsize,const string& unitig,bool stranded, int position1, int position2){
	string pathBegin(getPathBegin(numberBegin));
	//	pathBegin=pathBegin.substr(0,begsize);
	if(pathBegin.empty() and numberBegin.size()!=0){
		cout<<"fail to recompose path begin"<<endl;
		return "";
		//		exit(0);
	}
	string pathEnd(getPathEnd(numberEnd));
	//	pathEnd=pathEnd.substr(0,endsize);
	if(pathEnd.empty() and numberEnd.size()!=0){
		cout<<"fail to recompose path end"<<endl;
		return "";
		//		exit(0);
	}
	string path;
	if(!numberBegin.empty()){
		path=(compactionBegin(unitig,pathBegin,kgraph));
		if(path.empty()){
			//			path=(compactionBegin(reversecomplement(unitig),(pathBegin),kgraph));
			if(path.empty()){
				cout<<unitig.substr(0,kgraph)<<endl
				<<unitig.substr(unitig.size()-kgraph,kgraph)<<endl
				<<pathBegin.substr(0,kgraph)<<endl
				<<pathBegin.substr(pathBegin.size()-kgraph,kgraph)<<endl;
				cout<<reversecomplement(unitig.substr(0,kgraph))<<endl
				<<reversecomplement(unitig.substr(unitig.size()-kgraph,kgraph))	<<endl
				<<reversecomplement(pathBegin.substr(0,kgraph))<<endl
				<<reversecomplement(pathBegin.substr(pathBegin.size()-kgraph,kgraph))<<endl;
				cout<<"fail compactbegin in recover"<<endl;
				failedCompaction++;
				path=(compactionEnd(pathBegin,unitig,kgraph));
				cout<<"path"<<path<<endl;
				//				exit(0);
				return "";
			}
		}
	}else{
		path=unitig;
		//		if(stranded){
		//			if(position1>0){
		//				path=unitig;
		//			}else{
		//				if(unitig.size()+position1<kgraph){
		//					path=unitig.substr(unitig.size()-kgraph);
		//				}else{
		//					path=unitig.substr(-position1);
		//				}
		//			}
		//		}else{
		//			if(position2>0){
		//				path=unitig;
		//			}else{
		//				if(unitig.size()+position2<kgraph){
		//					path=unitig.substr(unitig.size()-kgraph);
		//				}else{
		//					path=unitig.substr(-position2);
		//				}
		//			}
		//		}
	}
	if(!numberEnd.empty()){
		string final(compactionEnd(path, pathEnd, kgraph));
		if(!final.empty()){
			path=final;
		}else{
			final=(compactionEnd(reversecomplement(path), pathEnd, kgraph));
			if(final.empty()){
				cout<<"fail compactend..."<<endl;
				//				cout<<path<<endl;
				//				cout<<(pathEnd)<<endl;
				//				cout<<path.substr(0,kgraph)<<endl
				//				<<path.substr(path.size()-kgraph,kgraph)<<endl
				//				<<pathEnd.substr(0,kgraph)<<endl
				//				<<pathEnd.substr(pathEnd.size()-kgraph,kgraph)<<endl;
				//				cout<<reversecomplement(path.substr(0,kgraph))<<endl
				//				<<reversecomplement(path.substr(path.size()-kgraph,kgraph))	<<endl
				//				<<reversecomplement(pathEnd.substr(0,kgraph))<<endl
				//				<<reversecomplement(pathEnd.substr(pathEnd.size()-kgraph,kgraph))<<endl;
				failedCompaction++;
				return "";
				//						exit(0);
			}
		}
	}else{
		//		path=path.substr(0,readsize);
	}
	if(path.empty()){
		cout<<"Cant recompose path...."<<endl;
	}
	return path;
}



void MappingSupervisor::MapFromUnitigsErrors(const string& unitig){
	//pas besoin de map here !
	unordered_map<rNumber,uint32_t> Candidate;
	unordered_map<rNumber,unordered_set<minimizer>> read2min;
	findCandidate(unitig,Candidate,read2min);
	if(Candidate.empty()){return;}
	bool done(false);
	unordered_multimap<string,string> genomicKmers(allKmerMapStranded(k2,unitig,nuc));
	vector<path> ListSons(listPathSons(offset, unitig.substr(unitig.size()-kgraph,kgraph),0));
	vector<path> ListFathers(listPathSons(offset, reversecomplement(unitig.substr(0,kgraph)),0));
	if(ListFathers.empty() and ListSons.empty()){
		island++;
		return;
	}
	string read,beg,end,path,seq1,seq2;
	vector<uNumber> numberBegin,numberEnd;
	//foreach reads that could map on the unitig (prepremapped)
	for(auto it=Candidate.begin(); it!=Candidate.end(); ++it){
		if(it->second>=multi){
			rNumber readnumber(it->first);
			candidateNumber++;
			if(number2position[readnumber]==0){
				continue;
			}
			read=getRead(readnumber);
			int position,positionUnitig,position1, position2;
			bool stranded;
			if(preMapUnitig(unitig, read, genomicKmers, position, read2min[it->first], positionUnitig, stranded,position1,position2)){
				if(!done){
					done=true;
					unitigsPreMapped++;
				}
				numberBegin=numberEnd={};
				beg=end="";
				if(mapOnGraph(numberBegin, numberEnd, unitig, read, position, ListFathers, ListSons, readnumber,beg,end)){
					path=(recoverPath(numberBegin, numberEnd, beg.size(), end.size(), read.size(), unitig, stranded, position1, position2));
					if(path.empty()){
						cout<<read<<endl;
						break;
					}
					seq1=seq2="";
					if(numberEnd.empty() and numberBegin.empty()){
						break;
						//						if(checking){
						//							positionUnitig=max(0,-positionUnitig);
						//							nw(unitig.substr(positionUnitig,read.size()), read, seq1, seq2, false);
						//						}
						//						if(!checking or scoreFromAlignment(seq1,seq2)<=errorRate){
						//							auto pathKmers(allKmerMapStranded(k2,unitig,nuc));
						//							double score(jaccardStrandedErrors(k2,read,pathKmers,nuc));
						//							globalscore=globalscore+score;
						//							if(readused.count(readPosition==0)){
						//								outFile<<">"<<++pathNumber<<endl<<path<<endl;
						//								readused.insert(readPosition);
						//								readInUnitig++;
						//								regionmapped+=read.size();
						//							}
						//						}else{
						//							fail++;
						//						}
					}else{
						if(checking){
							nw(path, read, seq1, seq2, false);
						}
						if(!checking or scoreFromAlignment(seq1,seq2)<=errorRate){
							//							auto pathKmers(allKmerMapStranded(k2,path,nuc));
							//							double score(jaccardStrandedErrors(k2,read,pathKmers,nuc));
							//							globalscore=globalscore+score;
							pathlength+=numberBegin.size()+numberEnd.size()+1;
							//							if(readused.count(readPosition)==0){
							//								readused.insert(readPosition);
							mutexEraseReads.lock();
							number2position[readnumber]=0;
							outFile<<">"<<++pathNumber<<endl<<path<<endl;
							mutexEraseReads.unlock();
							++aligneOnPathSucess;
							regionmapped+=read.size();

						}else{fail++;}
					}
					continue;
				}
			}
		}
	}
}






void MappingSupervisor::MapPart(){
	//foreach unitig (sort of)
	string unitig;
	while(true){
		mutexunitig.lock();
		if(indice<unitigs.size()){
			unitig=unitigs[indice++];
			mutexunitig.unlock();
		}else{
			mutexunitig.unlock();
			break;
		}
		if(unitig.size()>=minSizeUnitigs){
			++bigUnitig;
			if(bigUnitig%1000==0){
				cout<<bigUnitig<<endl;
				cout<<"Unitigs mapped "<<unitigsPreMapped<<" Percent unitigs mapped : "<<(100*unitigsPreMapped)/(bigUnitig)<<endl;
				cout<<"Read pre-mapped : "<<readMapped<<" Percent read pre-mapped : "<<(100*readMapped)/(readNumber)<<endl;
				cout<<"Read mapped on graph: "<<aligneOnPathSucess<<" Percent read mapped  on graph: "<<(100*aligneOnPathSucess)/(readNumber-readInUnitig+1)<<endl;
				cout<<"Read mapped : "<<aligneOnPathSucess+readInUnitig<<" Percent read mapped  "<<(100*aligneOnPathSucess+100*readInUnitig)/(readNumber+1)<<endl;
				cout<<island<<" islands..."<<endl;
				cout<<regionmapped/(1000*1000)<<" Mnuc mapped or "<<regionmapped/(1000*1)<<" Knuc "<<endl;
				cout<<leftmap<<" left map "<< rightmap<<" right map"<<endl;
				cout<<leftmapFail<<" left failt "<<rightmapFail<<" right fail"<<endl;
				cout<<readInUnitig<<" reads in unitig"<<endl;
				if(fail+aligneOnPathSucess!=0){
					cout<<"%sucess of the shared kmer heuristic : "<<(100*aligneOnPathSucess)/(fail+aligneOnPathSucess+1)<<endl;;
				}
				cout<<"mean score : "<<globalscore/(aligneOnPathSucess+1)<<endl;
				cout<<"fail from depth : "<<deepper<<endl;
				cout<<"failed compaction(bug) : "<<failedCompaction<<endl;
				cout<<"mean number of candidate: "<<candidateNumber/(bigUnitig+1)<<" "<<candidateNumber<<endl;
				cout<<"mean number of pathlent: "<<pathlength/(aligneOnPathSucess+1)<<endl;
				cout<<endl;
			}

			if(unitig.size()<offset){
				vector<path> ListFathers(listPathFathers(offset, unitig.substr(0,kgraph),0));
				for(size_t j(0); j<ListFathers.size(); ++j){
					path pathFather(ListFathers[j]);
					string newUnitig(compactionBegin(unitig, pathFather.str, kgraph));
					MapFromUnitigsErrors(newUnitig);
				}
				vector<path> ListSons(listPathSons(offset, unitig.substr(unitig.size()-kgraph,kgraph),0));
				for(size_t j(0); j<ListSons.size(); ++j){
					path pathSon(ListSons[j]);
					string newUnitig(compactionEnd(unitig, pathSon.str, kgraph));
					MapFromUnitigsErrors(newUnitig);
				}
			}else{
				MapFromUnitigsErrors(unitig);
			}
		}
	}
}



string MappingSupervisor::getPathEnd(vector<uNumber>& numbers){
	if(numbers.empty()){
		return "";
	}

	string path(unitigs[numbers[0]]);
	//	cout<<path<<endl;
	//	cout<<"getpathEndunitig"<<endl;
	//	cout<<path<<endl;
	//	cout<<numbers[0]<<" "<<endl;
	for(size_t i(1); i<numbers.size(); ++i){
		//				cout<<numbers[i]<<" "<<endl;
		string unitig(unitigs[numbers[i]]);
		//		cout<<unitig<<endl;
		string inter=compactionEnd(path, unitig, kgraph);
		if(inter.empty()){
			path=compactionEnd(reversecomplement(path), unitig, kgraph);
			if(path.empty()){
				cout<<"!"<<path<<"!"<<unitig<<endl;
//				cin.get();
			}
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
	//	reverse(numbers.begin(), numbers.end());
	string path(unitigs[numbers[0]]);

	for(size_t i(1); i<numbers.size(); ++i){
		string unitig(unitigs[numbers[i]]);
		string inter=compactionBegin(path, unitig, kgraph);
		if(inter.empty()){
			path=compactionBegin(reversecomplement(path), unitig, kgraph);
			if(path.empty()){
				cout<<"wtf : "<<path<<" "<<unitig<<endl;
				//				cin.get();
			}
		}else{
			path=inter;
		}
	}
	return path;
}





void MappingSupervisor::MapAll(){
	vector<thread> threads;

	for(size_t i(0); i<nbThreads; ++i) {
		threads.push_back(thread(&MappingSupervisor::MapPart, this));
	}

	for(auto &t : threads){t.join();}

	cout<<bigUnitig<<endl;
	cout<<"Unitigs mapped "<<unitigsPreMapped<<" Percent unitigs mapped : "<<(100*unitigsPreMapped)/(bigUnitig+1)<<endl;
	cout<<"Read pre-mapped : "<<readMapped<<" Percent read pre-mapped : "<<(100*readMapped)/(readNumber)<<endl;
	cout<<"Read mapped on graph: "<<aligneOnPathSucess<<" Percent read mapped  on graph: "<<(100*aligneOnPathSucess)/(readNumber-readInUnitig+1)<<endl;
	cout<<"Read mapped : "<<aligneOnPathSucess+readInUnitig<<" Percent read mapped  "<<(100*aligneOnPathSucess+100*readInUnitig)/(readNumber+1)<<endl;
	cout<<island<<" islands..."<<endl;
	cout<<regionmapped/(1000*1000)<<" Mnuc mapped or "<<regionmapped/(1000*1)<<" Knuc "<<endl;
	cout<<leftmap<<" left map "<< rightmap<<" right map"<<endl;
	cout<<leftmapFail<<" left failt "<<rightmapFail<<" right fail"<<endl;
	cout<<readInUnitig<<" reads in unitig"<<endl;
	if(fail+aligneOnPathSucess!=0){
		cout<<"%sucess of the shared kmer heuristic : "<<(100*aligneOnPathSucess)/(fail+aligneOnPathSucess+1)<<endl;;
	}
	cout<<"mean score : "<<globalscore/(aligneOnPathSucess+1)<<endl;
	cout<<"fail from depth : "<<deepper<<endl;
	cout<<"failed compaction(bug) : "<<failedCompaction<<endl;
	cout<<"mean number of candidate: "<<candidateNumber/(bigUnitig+1)<<" "<<candidateNumber<<endl;
	cout<<"mean number of pathlent: "<<pathlength/(aligneOnPathSucess+1)<<endl;
	cout<<endl;
}





