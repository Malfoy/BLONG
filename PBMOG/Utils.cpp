//
//  Utils.cpp
//  PBMOG
//
//  Created by malfoy on 09/01/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#include "Utils.h"

using namespace std;


string reversecompletment(const string& str){
	string res(str);
	int n = (int)str.size();
	for(int i(n-1), j(0); i > -1; i--, j++)
	{
		unsigned char c = str[i];
		unsigned char d = (c >> 4)&7;
		if (d >= 6) //(c >= 'a' && c <= 't')
		{
			// translates acgt to tgca
			c ^= 4;
			if ((c&3) != 3)
				c ^= 17;
			res[j] = c;
			continue;
		}
		if (d == 2)  // switch '+' with '-'
			res[j] = c ^ 6;
		else
		{
			// else it will only be a number, just copy it
			res[j] = c;
		}
	}
	return res;
}

char nuc2int(char c){
	switch(c){
		case 'a': return 0;
		case 'c': return 1;
		case 'g': return 2;
		case 't': return 3;
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
	}
	cout<<"fuck"<<endl;
	return 0;
}

uint64_t xorshift64(uint64_t x) {
	x ^= x >> 12; // a
	x ^= x << 25; // b
	x ^= x >> 27; // c
	return x * UINT64_C(2685821657736338717);
}

vector<minimizer> allHash(size_t k,const string& seq){
	vector<minimizer> sketch;
	for(size_t i(0);i+k<seq.size();++i){
		sketch.push_back(seq2int(seq.substr(i,k)));
	}
	return sketch;
}

double jaccard3(size_t k, const string& seq,const unordered_set<minimizer>& A){
	minimizer kmer;
	double inter(0);

	for(size_t i(0);i+k<seq.size();++i){
		kmer=seq2int(seq.substr(i,k));
		if(A.count(kmer)>0){
			++inter;
		}
	}
	return double(100*inter/(A.size()));
}

void minHash2(size_t H, size_t k, const string& seq, vector<minimizer>& previous){
	vector<uint64_t> sketch(H);
	vector<minimizer> sketchs(H);
	uint64_t hashValue;
	uint32_t kmer;
	//~ hash<uint32_t> hash;

	kmer=seq2int(seq.substr(0,k));
	//~ hashValue=hash(kmer);
	hashValue=xorshift64(kmer);
	for(size_t j(0);j<H;++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=xorshift64(hashValue);
	}
	for(size_t i(1);i+k<seq.size();++i){
		kmer=seq2int(seq.substr(i,k));
		hashValue=xorshift64(kmer);
		//~ hashValue=hash(kmer);
		for(size_t j(0);j<H;++j){
			if(hashValue<sketch[j]){
				sketch[j]=hashValue;
				sketchs[j]=kmer;
			}
			hashValue=xorshift64(hashValue);
		}
	}
	previous.insert(previous.end(),sketchs.begin(),sketchs.end());
}


uint32_t seq2int(const string& seq){
	uint32_t res(0);
	for(uint i(0);i<seq.size();++i){
		res+=nuc2int(seq[i]);
		res<<=2;
	}
	//~ cout<<seq<<" "<<res<<endl;
	string rc(reversecompletment(seq));
	uint32_t res2(0);
	for(uint i(0);i<rc.size();++i){
		res2+=nuc2int(rc[i]);
		res2<<=2;
	}

	return min(res,res2);
}

vector<minimizer> minHashpart(size_t H, size_t k,const string& seq, size_t part){
	vector<minimizer> result;
	size_t size(seq.size()/part);
	for(size_t i(0);i<part;++i){
		minHash2(H/part,k,seq.substr(i*size,size+k),result);
	}
	return result;
}
