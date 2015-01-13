//
//  Utils.h
//  PBMOG
//
//  Created by malfoy on 09/01/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#ifndef __PBMOG__Utils__
#define __PBMOG__Utils__

#include <stdio.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <functional>
#include <unordered_set>
//#include <sparsehash/sparse_hash_map>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <atomic>



//~ #define unordered_map sparse_hash_map
#define minimizer uint32_t
#define rNumber uint32_t

using namespace std;

uint32_t seq2int(const string& seq);
string reversecompletment(const string& str);
char nuc2int(char c);
vector<minimizer> minHashpart(size_t H, size_t k,const string& seq, size_t part);
void minHash2(size_t H, size_t k, const string& seq, vector<minimizer>& previous);
uint64_t xorshift64(uint64_t x);
vector<minimizer> allHash(size_t k,const string& seq);
double jaccard3(size_t k, const string& seq,const unordered_set<minimizer>& A);

#endif /* defined(__PBMOG__Utils__) */
