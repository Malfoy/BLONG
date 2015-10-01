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
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <unordered_map>
#include <functional>
#include <unordered_set>
#include <set>
#include <sparsehash/dense_hash_map>
#include <sparsehash/dense_hash_set>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <cctype>
#include <atomic>
#include <random>


//~ #define unordered_map sparse_hash_map
//#define unordered_map dense_hash_map
#define minimizer uint32_t
#define rNumber uint32_t
#define uNumber uint32_t


using namespace std;
using namespace google;

minimizer seq2int(const string& seq);
minimizer seq2intStranded(const string& seq);
string reversecomplement(const string& str);
char nuc2int(char c);
vector<minimizer> minHashpart(size_t H, size_t k,const string& seq, size_t part);
void minHash2(size_t H, size_t k, const string& seq, vector<minimizer>& previous);
uint64_t xorshift64(uint64_t x);
vector<minimizer> allHash(size_t k,const string& seq);
double jaccard(size_t k, const string& seq,const unordered_set<minimizer>& A);
double jaccardStranded(size_t k, const string& seq,const unordered_set<minimizer>& A);
double jaccardAlt(size_t k, const string& seq,const unordered_set<minimizer>& A);
vector<size_t> bounds(size_t n,size_t size);
string getRepresent(const string& str);
string compaction(const string& seq1,const  string& seq2, size_t k);
string compactionEnd(const string& seq1,const  string& seq2, size_t k);
string compactionBegin(const string& seq1,const  string& seq2, size_t k);
void readContigsforstats(const string& File, size_t k, bool elag, bool compact,bool unitigb);
unordered_map<string,vector<minimizer>> getGraph(const vector<string>& unitigs, size_t k);
vector<string> loadFASTQ(const string& unitigFile,bool homo,size_t size,char frac);
vector<string> loadFASTA(const string& unitigFile,bool homo,size_t size, size_t frac);
string homocompression(const string& seq);
vector<string> loadUnitigs(const string& unitigFile,bool homo);
void minHash3(size_t H, size_t k,const string& seq, vector<minimizer>& previous, const unordered_set<minimizer>& filter);
uint64_t xorshift(uint64_t x);
vector<minimizer> minHashpart2(size_t H, size_t k,const string& seq, size_t part, const unordered_set<minimizer>& filter);
unordered_set <minimizer> allKmerSet(size_t k,const string& seq);
unordered_set <minimizer> allKmerSetStranded(size_t k,const string& seq);
int positionInSeq(const string& seq, minimizer min, size_t k);
int positionInSeqStranded(const string& seq, minimizer min, size_t k);
int positionInSeqStrandedEnd(const string& seq, minimizer min, size_t k);
size_t random(size_t max);
string randomString( size_t length );
void updateMinimizer(minimizer& min, char nuc,size_t k);
void updateMinimizerRC(minimizer&	min, char nuc,size_t k);
void updateMinimizerEnd(minimizer&	min, char nuc,size_t k);
vector<string> kmerCounting(const string& fileName,size_t k);
unordered_multimap<string,string> allKmerMapStranded(size_t k,const string& seq,char nuc);
double jaccardStrandedErrors(size_t k, const string& seq, const unordered_multimap<string, string>& genomicKmers,char nuc);
void printMinimizer(minimizer min,size_t k);
bool isCorrect(const string& seq,const string& ref);
double scoreFromAlignment(const string& seq1,const string& seq2);


#endif /* defined(__PBMOG__Utils__) */
