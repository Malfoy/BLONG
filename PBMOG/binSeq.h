//
//  binSeq.h
//  PBMOG
//
//  Created by malfoy on 27/01/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#ifndef __PBMOG__binSeq__
#define __PBMOG__binSeq__

#include <stdio.h>
#include <vector>

using namespace std;


class binSeq {
public:
	static const unsigned char rc[];
	vector<unsigned char> vect;

	string str();//
	binSeq sub(size_t begin);
	uint64_t getBegin(size_t size);
	uint64_t getEnd(size_t size);
	void reverse();
	void add(binSeq);
	void resize();

	binSeq(const string& str);
	binSeq();
	binSeq(const binSeq& bs);

};



char char2int(char c);
char int2char(char c);


#endif /* defined(__PBMOG__binSeq__) */
