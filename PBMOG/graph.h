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
	unordered_map <string,vector<uNumber>> gotRight;
	unordered_map <string,vector<uNumber>> gotLeft;

	void addRight(const string& str, uNumber i){
		gotRight[str].push_back(i);
	}
	void addLeft(const string& str, uNumber i){
		gotLeft[str].push_back(i);

	}

	vector<uNumber> getRight(string str){
		str=getRepresent(str);
		if(gotRight.count(str)!=0){
			return gotRight[str];
		}
		return {};
	}

	vector<uNumber> getLeft(string str){
		str=getRepresent(str);
		if(gotLeft.count(str)!=0){
			return gotLeft[str];
		}
		return {};
	}

	graph(size_t kgraph){
		gotRight.set_empty_key("0");
		gotLeft.set_empty_key("0");
		k=kgraph;
	};

	graph(){
	};

	graph(const vector<string>& unitigs,size_t kgraph){
		gotRight.set_empty_key("0");
		gotLeft.set_empty_key("0");
		k=kgraph;
		string unitig,seq1,seq2;
		for(uNumber i(0);i<unitigs.size();++i){
			unitig=(unitigs[i]);
			seq1=getRepresent(unitig.substr(0,k));
			seq2=getRepresent(unitig.substr(unitig.size()-k,k));
			addLeft(seq1, i);
			addRight(seq2, i);
		}
	};

	}

	;

#endif /* defined(__PBMOG__graph__) */
