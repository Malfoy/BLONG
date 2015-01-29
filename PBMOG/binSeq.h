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
	bool isNumber;

	string str();
	binSeq sub(size_t begin);
	binSeq sub(size_t begin,size_t size);
	uint64_t getBegin(size_t size);
	uint64_t getEnd(size_t size);
	void reverse();
	void add(binSeq);
	void resize();
	size_t size();
	uint32_t getInt();

	binSeq& operator+=(const binSeq& rhs){
		this->add(rhs);
		return *this;
	}


	

	binSeq(const string& str);
	binSeq();
	binSeq(const binSeq& bs);
	binSeq(uint32_t);

};

inline binSeq operator+(binSeq lhs, const binSeq& rhs){
	lhs += rhs;
	return lhs;
}



unsigned char char2int(unsigned char c);
unsigned char int2char(unsigned char c);

void testBinSeq();


#endif /* defined(__PBMOG__binSeq__) */
