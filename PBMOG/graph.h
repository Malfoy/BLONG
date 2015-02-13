//
//  graph.h
//  PBMOG
//
//  Created by malfoy on 13/02/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#ifndef __PBMOG__graph__
#define __PBMOG__graph__

#include <stdio.h>
#include "Utils.h"


class graph{
public:
	size_t k;
	unordered_map <string,vector<uNumber>> right;
	unordered_map <string,vector<uNumber>> left;

	void addBeg(const string& str, uNumber i){
		string rc(reversecomplement(str));
		if(str<rc){
			left[str].push_back(i);
		}else{
			right[rc].push_back(i);
		}
	}

	void addEnd(const string& str, uNumber i){
		string rc(reversecomplement(str));
		if(str<reversecomplement(str)){
			right[str].push_back(i);
		}else{
			left[rc].push_back(i);
		}
	}


	vector<uNumber> getBegin(string str){
		string rc(reversecomplement(str));
		if(str<rc){
			if(left.count(str)!=0){
				return left[str];
			}else{
				return {};
			}
		}else{
			if(right.count(rc)!=0){
				return right[rc];
			}else{
				return {};
			}
		}
	}


	vector<uNumber> getEnd(string str){
		string rc(reversecomplement(str));
		if(str<rc){
			if(right.count(str)!=0){
				return right[str];
			}else{
				return {};
			}
		}else{
			if(left.count(rc)!=0){
				return left[rc];
			}else{
				return {};
			}
		}
	}

//	vector<uNumber> getLeft(string str){
//		str=getRepresent(str);
//		if(gotLeft.count(str)!=0){
//			return gotLeft[str];
//		}
//		return {};
//	}

	graph(size_t kgraph){
		right.set_empty_key("0");
		left.set_empty_key("0");
		k=kgraph;
	};

	graph(){
	};

	graph(const vector<string>& unitigs,size_t kgraph){
		right.set_empty_key("0");
		left.set_empty_key("0");
		k=kgraph;
		string unitig,seq1,seq2;
		for(uNumber i(0);i<unitigs.size();++i){
			unitig=(unitigs[i]);
			seq1=(unitig.substr(0,k));
			seq2=(unitig.substr(unitig.size()-k,k));
			addBeg(seq1,i);
			addEnd(seq2,i);
		}
	};

	}

	;

#endif /* defined(__PBMOG__graph__) */
