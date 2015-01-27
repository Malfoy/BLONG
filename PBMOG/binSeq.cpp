//
//  binSeq.cpp
//  PBMOG
//
//  Created by malfoy on 27/01/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#include "binSeq.h"
#include <string>



const unsigned char binSeq::rc[]={
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,

	0b01000011,0b01000010,0b01000001,0b01000000,0b01000000,0b01000000,0b01000000,0b01000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,

	0b10001111,0b10001011,0b10000111,0b10000011,0b10001110,0b10001010,0b10000110,0b10000010,
	0b00001101,0b00001001,0b00000101,0b00000001,0b00001100,0b00001000,0b00000100,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,

	0b11111111,0b11101111,0b11011111,0b11001111,0b11111011,0b11101011,0b11011011,0b11001011,
	0b11110111,0b11100111,0b11010111,0b11000111,0b11110011,0b11100011,0b11010011,0b11000011,
	0b11111110,0b11101110,0b11011110,0b11001110,0b11111010,0b11101010,0b11011010,0b11001010,
	0b11110110,0b11100110,0b11010110,0b11000110,0b11110010,0b11100010,0b11010010,0b11000010,
	0b11111101,0b11101101,0b11011101,0b11001101,0b11111001,0b11101001,0b11011001,0b11001001,
	0b11110101,0b11100101,0b11010101,0b11000101,0b11110001,0b11100001,0b11010001,0b11000001,
	0b11111100,0b11101100,0b11011100,0b11001100,0b11111000,0b11101000,0b11011000,0b11001000,
	0b11110100,0b11100100,0b11010100,0b11000100,0b11110000,0b11100000,0b11010000,0b11000000,

};



char char2int(char c){
	switch (c) {
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
	}
	return 3;
}



char int2char(char c){
	switch (c) {
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
	}
	return 'T';
}




binSeq::binSeq(const string& str){
	char mod(str.size()%3);
	char c;
	for (size_t i(0); i<str.size()-mod; i+=3){
		c=12;
		c+=char2int(str[i]);
		c<<=2;
		c+=char2int(str[i+1]);
		c<<=2;
		c+=char2int(str[i+2]);
		vect.push_back(c);
	}

	switch (mod) {
		case 0:
			break;
		case 1:
			c=1<<6;
			c+=char2int(str[str.size()-1]);
			vect.push_back(c);
			break;
		case 2:
			c=2<<6;
			c+=char2int(str[str.size()-2]);
			c<<=2;
			c+=char2int(str[str.size()-1]);
			vect.push_back(c);
			break;
	}
}



binSeq::binSeq(const binSeq& bs){
	vect=bs.vect;
}



string binSeq::str(){
	string res;
	for(size_t i(0);i<vect.size();++i){
		char c(vect[i]);
		char mod(c%(1<<6));
		switch (mod) {
			case 1:
				res.push_back(int2char(c/4));
				break;
			case 2:
				res.push_back(int2char(c/4));
				c>>=2;
				res.push_back(int2char(c/4));
				break;

			case 3:
				res.push_back(int2char(c/4));
				c>>=2;
				res.push_back(int2char(c/4));
				c>>=2;
				res.push_back(int2char(c/4));
				break;
		}
	}
	return res;
}


binSeq binSeq::sub(size_t begin){
	binSeq res;
	size_t count(0);
	bool go(true);
	size_t i(0);
	for(; i<vect.size() and go; ++i){
		char c(vect[i]);
		char n(c%(1<<6));
		if(count+n<begin){
			count+=n;
		}else{
			if(count+n==begin){
				go=false;
			}else{
				go=false;
				char toGet((char)(count+n-begin)	);
				char nc((char)(toGet<<6));
				nc+=c%(1<<(2*toGet));
				res.vect.push_back(nc);
			}
		}
	}
	res.vect.insert(res.vect.end(), vect.begin()+i, vect.end());

	return res;
}



uint64_t binSeq::getBegin(size_t size){
	uint64_t res(0);
	size_t i(0);
	for(size_t c(0); c<size; ++c){
		res+=vect[i]%(1<<6);
		i+=vect[i]/(1<<6);
	}
	return res;
}



uint64_t binSeq::getEnd(size_t size){
	uint64_t res(0);
	size_t i(vect.size()-1);
	for(size_t c(0); c<size; ++c){
		res+=vect[i]%(1<<6);
		i-=vect[i]/(1<<6);
	}
	return res;
}



void binSeq::add(binSeq bs){
	vect.insert(vect.begin(), bs.vect.begin(), bs.vect.end());
}


void binSeq::reverse(){
	vector<unsigned char> V;
	for(int i((int)vect.size()-1); i>-1; --i){
		V.push_back(rc[vect[i]]);
	}
	vect=V;
}









































