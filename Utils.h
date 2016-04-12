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
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <cctype>
#include <atomic>
#include <random>
#include <mutex>


//~ #define unordered_map sparse_hash_map
//#define unordered_map dense_hash_map
#define minimizer uint32_t
#define rNumber uint32_t
#define rPosition uint64_t
#define uNumber uint32_t


using namespace std;

int countFASTA(const string& seqFile,uint  sizeMin);
minimizer seq2int(const string& seq);
minimizer seq2intStranded(const string& seq);
string reversecomplement(const string& str);
char nuc2int(char c);
vector<minimizer> minHashpart(uint  H, uint  k,const string& seq, uint  part);
void minHash2(uint  H, uint  k, const string& seq, vector<minimizer>& previous);
uint64_t xorshift64(uint64_t x);
vector<minimizer> allHash(uint  k,const string& seq);
double jaccard(uint  k, const string& seq,const unordered_set<minimizer>& A);
double jaccardStranded(uint  k, const string& seq,const unordered_set<minimizer>& A);
double jaccardAlt(uint  k, const string& seq,const unordered_set<minimizer>& A);
vector<uint > bounds(uint  n,uint  size);
string getRepresent(const string& str);
string compaction(const string& seq1,const  string& seq2, uint  k);
string compactionEnd(const string& seq1,const  string& seq2, uint  k);
string compactionBegin(const string& seq1,const  string& seq2, uint  k);
void readContigsforstats(const string& File, uint  k, bool elag, bool compact,bool unitigb);
unordered_map<string,vector<minimizer>> getGraph(const vector<string>& unitigs, uint  k);
vector<string> loadFASTQ(const string& unitigFile,bool homo,uint  size,char frac);
vector<string> loadFASTA(const string& unitigFile,bool homo,uint  size, uint  frac);
string homocompression(const string& seq);
vector<string> loadUnitigs(const string& unitigFile,bool homo);
void minHash3(uint  H, uint  k,const string& seq, vector<minimizer>& previous, const unordered_set<minimizer>& filter);
uint64_t xorshift(uint64_t x);
vector<minimizer> minHashpart2(uint  H, uint  k,const string& seq, uint  part, const unordered_set<minimizer>& filter);
unordered_set <minimizer> allKmerSet(uint  k,const string& seq);
unordered_set <minimizer> allKmerSetStranded(uint  k,const string& seq);
int positionInSeq(const string& seq, minimizer min, uint  k);
int positionInSeqStranded(const string& seq, minimizer min, uint  k);
int positionInSeqStrandedEnd(const string& seq, minimizer min, uint  k);
uint  random(uint  max);
string randomString( uint  length );
void updateMinimizer(minimizer& min, char nuc,uint  k);
void updateMinimizerRC(minimizer&	min, char nuc,uint  k);
void updateMinimizerEnd(minimizer&	min, char nuc,uint  k);
vector<string> kmerCounting(const string& fileName,uint  k);
unordered_multimap<uint32_t,uint32_t> allKmerMapStranded(const char k,const string& seq, const  char nuc);
double jaccardStrandedErrors(uint  k, const string& seq, const unordered_multimap<string, string>& genomicKmers,char nuc);
void printMinimizer(minimizer min,uint  k);
bool isCorrect(const string& seq,const string& ref);
double scoreFromAlignment(const string& seq1,const string& seq2);
void removeDuplicate(vector<minimizer>& vec);
unordered_map<minimizer,vector<rNumber>> indexSeqDisk(const string& seqs, uint  H, uint  k, uint  part,vector<rPosition>& number2position,unordered_set<minimizer>& lol);
unordered_set<minimizer> allKmersetu(uint  k,const string& seq);
uint32_t jaccardStrandedErrors(char k, const string& seq, const unordered_multimap<uint32_t, uint32_t>& genomicKmers, char nuc);
void printMinimizer(minimizer min,uint  k);
vector<uint32_t>* allKmerVectStranded(const char k,const string& seq, const  char nuc);
uint32_t jaccardStrandedErrors(char k, const string& seq, vector<uint32_t>* genomicKmers, char nuc);
void updateMinimizer16(minimizer&	min, char nuc,uint  k);
unordered_set<minimizer> filterUnitigs(const vector<string>& unitigs, uint  k, uint  H, uint  part);



#endif /* defined(__PBMOG__Utils__) */
