#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <functional>
#include <unordered_set>
#include <sparsehash/sparse_hash_map>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>

//~ #define unordered_map sparse_hash_map
#define minimizer uint32_t
#define rNumber uint32_t

using namespace std;
using namespace google;


mutex myMutex;



//~ uint64_t xorshift(uint64_t x) {
//~ x ^= x >> 12; // a
//~ x ^= x << 25; // b
//~ x ^= x >> 27; // c
//~ return x * UINT64_C(2685821657736338717);
//~ }


string bwt(string input,size_t &primaryidx) {
 string bwt_str;
	vector<string> rotations;

	// 1. generate rotations of input
	for(uint n=0;n<input.size();n++) {
		string s;

		s += input[n];
		for(uint i=n+1;i!=n;i++) {
			if(i==input.size()) i=0;
			s += input[i];
			if(i==input.size()-1) i=-1;
		}
		rotations.push_back(s);
	}

	// 2. sort
	sort(rotations.begin(),rotations.end());

	bwt_str.clear();
	for(size_t n=0;n<rotations.size();n++) {
		bwt_str += rotations[n].substr(rotations[n].size()-1,1);
		if(rotations[n] == input) {primaryidx = n;}
	}
	return bwt_str;
}




string rbwt(string bwt_str,int primaryidx) {

	vector<string> cur_str;

	for(size_t n=0;n<bwt_str.size();n++) {
		string s;
		s = bwt_str.substr(n,1);
		cur_str.push_back(s);
	}

	for(;cur_str[0].size() < bwt_str.size();) {
		vector<string> new_str = cur_str;
		sort(new_str.begin(),new_str.end());

		for(size_t n=0;n<cur_str.size();n++) {
			cur_str[n] = cur_str[n] + new_str[n].substr(new_str[n].size()-1,1);
		}
	}

	//~ for(size_t n=0;n<cur_str.size();n++) {
	//~ cout << cur_str[n] << endl;
	//~ }

	sort(cur_str.begin(),cur_str.end());


	return cur_str[primaryidx];

}




string homocompression(const string& seqo){
	string res;
	size_t ind;
	string seq=bwt(seqo,ind);
	res.push_back(seq[0]);
	char last(seq[0]);
	for(uint i(1);i<seq.size();++i){
		if(seq[i]!=last){
			res.push_back(seq[i]);
			last=seq[i];
		}
	}
	//~ cout<<seq<<endl<<res<<endl;
	//~ cin.get();
	return res;
}




int myrandom (int i) { return rand()%i;}

//static void PrintAlignment(const StripedSmithWaterman::Alignment& alignment){
//	cout << "===== SSW result =====" << endl;
//	cout << "Best Smith-Waterman score:\t" << alignment.sw_score << endl
//	<< "Next-best Smith-Waterman score:\t" << alignment.sw_score_next_best << endl
//	<< "Reference start:\t" << alignment.ref_begin << endl
//	<< "Reference end:\t" << alignment.ref_end << endl
//	<< "Query start:\t" << alignment.query_begin << endl
//	<< "Query end:\t" << alignment.query_end << endl
//	<< "Next-best reference end:\t" << alignment.ref_end_next_best << endl
//	<< "Number of mismatches:\t" << alignment.mismatches << endl
//	<< "Cigar: " << alignment.cigar_string << endl;
//	cout << "======================" << endl;
//	cin.get();
//}


//
//bool score(const string& query, const string& ref){
//	//~ return true;
//
//	StripedSmithWaterman::Aligner aligner;
//	StripedSmithWaterman::Filter filter;
//	StripedSmithWaterman::Alignment alignment;
//	if(aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment)){
//		//~ if((alignment.query_end-alignment.query_begin)>0.95*query.size()){
//		//~ cout<<(alignment.query_end-alignment.query_begin)<<endl;
//		//~ cout<<query.size()<<endl;
//		//~ cin.get();
//		//~ cout<<(100*alignment.sw_score)/(query.size())<<endl;
//		if((100*alignment.sw_score)/(query.size())>85){
//			return true;
//		}else{
//			return false;
//		}
//		//~ }
//	}
//	return false;
//}

uint32_t xorshift(uint32_t x,uint32_t y=12,uint32_t z=25,uint32_t w=27){
	uint32_t t = x ^ (x << 11);
	x = y; y = z; z = w;
	return w = w ^ (w >> 19) ^ t ^ (t >> 8);
}




uint64_t xorshift64(uint64_t x) {
	x ^= x >> 12; // a
	x ^= x << 25; // b
	x ^= x >> 27; // c
	return x * UINT64_C(2685821657736338717);
}




string reversecompletment(const string& str){
	string res(str);
	size_t n = str.size();
	for(size_t i(n-1), j(0); i > -1; i--, j++)
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




string getRepresent(const string& str){
	return(min(str,reversecompletment(str)));
}




uint32_t murmur(const char *key, uint32_t len=16, uint32_t seed=12) {
	static const uint32_t c1 = 0xcc9e2d51;
	static const uint32_t c2 = 0x1b873593;
	static const uint32_t r1 = 15;
	static const uint32_t r2 = 13;
	static const uint32_t m = 5;
	static const uint32_t n = 0xe6546b64;

	uint32_t hash = seed;

	const int nblocks = len / 4;
	const uint32_t *blocks = (const uint32_t *) key;
	int i;
	for (i = 0; i < nblocks; i++) {
		uint32_t k = blocks[i];
		k *= c1;
		k = (k << r1) | (k >> (32 - r1));
		k *= c2;

		hash ^= k;
		hash = ((hash << r2) | (hash >> (32 - r2))) * m + n;
	}

	const uint8_t *tail = (const uint8_t *) (key + nblocks * 4);
	uint32_t k1 = 0;

	switch (len & 3) {
		case 3:
			k1 ^= tail[2] << 16;
		case 2:
			k1 ^= tail[1] << 8;
		case 1:
			k1 ^= tail[0];

			k1 *= c1;
			k1 = (k1 << r1) | (k1 >> (32 - r1));
			k1 *= c2;
			hash ^= k1;
	}

	hash ^= len;
	hash ^= (hash >> 16);
	hash *= 0x85ebca6b;
	hash ^= (hash >> 13);
	hash *= 0xc2b2ae35;
	hash ^= (hash >> 16);

	return hash;
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




uint64_t seq2int2(const string& seq){
	uint64_t res(0);
	for(uint i(0);i<seq.size();++i){
		res+=nuc2int(seq[i]);
		res<<=2;
	}
	//~ cout<<seq<<" "<<res<<endl;
	string rc(reversecompletment(seq));
	uint64_t res2(0);
	for(uint i(0);i<rc.size();++i){
		res2+=nuc2int(rc[i]);
		res2<<=2;
	}

	return min(res,res2);
}




vector<minimizer> allHash(size_t k,const string& seq){
	vector<minimizer> sketch;
	for(uint i(0);i+k<seq.size();++i){
		sketch.push_back(seq2int(seq.substr(i,k)));
	}
	return sketch;
}




vector<string> allHashS(int k,const string& seq){
	vector<string> sketch;
	for(uint i(0);i+k<seq.size();++i){
		sketch.push_back(getRepresent(seq.substr(i,k)));
	}
	return sketch;
}




vector<minimizer> minHash(int H, int k, const string& seq){
	vector<minimizer> sketchs(H);
	vector<uint32_t> sketch(H);
	uint32_t hashValue;
	minimizer kmer;

	kmer=seq2int(seq.substr(0,k));
	hashValue=xorshift(kmer);
	for(int j(0);j<H;++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=xorshift(hashValue);
	}

	for(unsigned int i(1);i+k<seq.size();++i){
		kmer=seq2int(seq.substr(i,k));
		hashValue=xorshift(kmer);
		for(int j(0);j<H;++j){
			if(hashValue<sketch[j]){
				sketch[j]=hashValue;
				sketchs[j]=kmer;
			}
			hashValue=xorshift(hashValue);
		}
	}
	return sketchs;
}




vector<string> minHashS(int H, int k, const string& seq){
	vector<string> sketchs(H);
	vector<uint64_t> sketch(H);
	uint64_t hashValue;
	string kmer;

	kmer=getRepresent(seq.substr(0,k));
	hashValue=xorshift64(seq2int2(kmer));
	for(int j(0);j<H;++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=xorshift64(hashValue);
	}

	for(unsigned int i(1);i+k<seq.size();++i){
		kmer=getRepresent(seq.substr(i,k));
		hashValue=xorshift64(seq2int2(kmer));
		for(int j(0);j<H;++j){
			if(hashValue<sketch[j]){
				sketch[j]=hashValue;
				sketchs[j]=kmer;
			}
			hashValue=xorshift64(hashValue);
		}
	}
	return sketchs;
}




vector<string> minHashSf(int H, int k, const string& seq, const unordered_set<string>& filter){
	vector<string> sketchs(H);
	vector<uint64_t> sketch(H);
	uint64_t hashValue;
	string kmer;

	kmer=getRepresent(seq.substr(0,k));
	hashValue=xorshift64(seq2int2(kmer));
	for(int j(0);j<H;++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=xorshift64(hashValue);
	}

	for(unsigned int i(1);i+k<seq.size();++i){
		kmer=getRepresent(seq.substr(i,k));
		if(filter.count(kmer)!=0){
			hashValue=xorshift64(seq2int2(kmer));
			for(int j(0);j<H;++j){
				if(hashValue<sketch[j]){
					sketch[j]=hashValue;
					sketchs[j]=kmer;
				}
				hashValue=xorshift64(hashValue);
			}
		}
	}
	return sketchs;
}




void minHash2(size_t H, size_t k, const string& seq, vector<minimizer>& previous){
	vector<uint32_t> sketch(H);
	vector<minimizer> sketchs(H);
	uint32_t hashValue;
	uint32_t kmer;
	//~ hash<uint32_t> hash;

	kmer=seq2int(seq.substr(0,k));
	//~ hashValue=hash(kmer);
	hashValue=xorshift(kmer);
	for(int j(0);j<H;++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=xorshift(hashValue);
	}

	for(unsigned int i(1);i+k<seq.size();++i){
		kmer=seq2int(seq.substr(i,k));
		hashValue=xorshift(kmer);
		//~ hashValue=hash(kmer);
		for(int j(0);j<H;++j){
			if(hashValue<sketch[j]){
				sketch[j]=hashValue;
				sketchs[j]=kmer;
			}
			hashValue=xorshift(hashValue);
		}
	}

	previous.insert(previous.end(),sketchs.begin(),sketchs.end());
}




void minHash3(size_t H, size_t k,const string& seq, vector<minimizer>& previous, const unordered_set<minimizer>& filter){
	vector<uint32_t> sketch(H);
	vector<minimizer> sketchs(H);
	uint32_t hashValue;
	uint32_t kmer;
	//~ hash<uint32_t> hash;

	kmer=seq2int(seq.substr(0,k));
	hashValue=xorshift(kmer);
	//~ hashValue=hash(kmer);
	for(int j(0);j<H;++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=xorshift(hashValue);
	}

	for(unsigned int i(1);i+k<seq.size();++i){
		kmer=seq2int(seq.substr(i,k));
		if(filter.count(kmer)!=0){
			hashValue=xorshift(kmer);
			//~ hashValue=hash(kmer);
			for(int j(0);j<H;++j){
				if(hashValue<sketch[j]){
					sketch[j]=hashValue;
					sketchs[j]=kmer;
				}
				hashValue=xorshift(hashValue);
			}
		}
	}
	previous.insert(previous.end(),sketchs.begin(),sketchs.end());
}




vector<minimizer> minHashpart(size_t H, size_t k,const string& seq, size_t part){
	vector<minimizer> result;
	size_t size(seq.size()/part);
	for(int i(0);i<part;++i){
		minHash2(H/part,k,seq.substr(i*size,size+k),result);
	}
	return result;
}




vector<minimizer> minHashpart2(size_t H, size_t k,const string& seq, size_t part, const unordered_set<minimizer>& filter){
	vector<minimizer> result;
	//~ cout<<"debut"<<endl;
	size_t size(seq.size()/part);
	for(size_t i(0);i<part;++i){
		//~ minHash3(H/part,k,seq.substr(i*size,size+16),result,filter);
		minHash3(H/part,k,seq.substr(i*size,size),result,filter);
	}
	//~ cout<<"lol"<<endl;
	return result;
}




int commonMin(int H, int k,const string& seq1,const string& seq2){
	int result(0);
	vector<uint32_t> sketch1=minHash(H,k,seq1);
	vector<uint32_t> sketch2=minHash(H,k,seq2);
	for(unsigned int i(0);i<sketch1.size();++i){
		if(sketch1[i]==sketch2[i]){
			++result;
		}
	}
	return result;
}




bool jaccard(int k, const string& seq1, const string& seq2, double score){
	unordered_set<uint32_t> A,B;
	uint32_t kmer;
	double inter(0),un(0);

	for(unsigned int i(0);i+k<seq1.size();++i){
		kmer=seq2int(seq1.substr(i,k));
		A.insert(kmer);

	}
	for(unsigned int i(0);i+k<seq2.size();++i){
		kmer=seq2int(seq2.substr(i,k));
		B.insert(kmer);

	}
	for (const uint32_t& x: A){
		++un;
		if(B.count(x)>0){
			++inter;
		}
	}

	//~ for (const int& x: B){
	//~ if(A.count(x)==0){
	//~ ++un;
	//~ }
	//~ }
	return (100*inter/un>score);
	//~ return inter;
}




double jaccard2(int k, const string& seq1, const string& seq2){
	unordered_set<minimizer> A,B;
	double inter(0),un(0);
	uint H(1000);
	vector<minimizer> sketch1,sketch2;

	//~ sketch1=allHash(k,seq1);
	sketch1=minHashpart(H,k,seq1,1);
	for(unsigned int i(0);i<sketch1.size();++i){
		A.insert(sketch1[i]);
	}


	sketch2=minHashpart(H,k,seq2,1);

	for(unsigned int i(0);i<sketch2.size();++i){
		B.insert(sketch2[i]);
	}

	for (const uint32_t& x: A){
		++un;
		if(B.count(x)>0){
			++inter;
		}
	}
	//~
	//~ for (const uint32_t& x: B){
	//~ ++un;
	//~ if(A.count(x)>0){
	//~ ++inter;
	//~ }
	//~ }

	//~ for(unsigned int i(0);i<sketch1.size();++i){
	//~ if(sketch1[i]==sketch2[i]){
	//~ inter++;
	//~ }
	//~ }

	//~ for (const uint32_t& x: B){
	//~ if(A.count(x)==0){
	//~ ++un;
	//~ }
	//~ }
	//~ cout<<inter<<" "<<un<<endl;
	return inter/un;
	//~ return inter;
}




double jaccard3(size_t k, const string& seq,const unordered_set<minimizer>& A){
	minimizer kmer;
	double inter(0);

	for(unsigned int i(0);i+k<seq.size();++i){
		kmer=seq2int(seq.substr(i,k));
		if(A.count(kmer)>0){
			++inter;
		}
	}
	return double(100*inter/(A.size()));
}




double jaccard4(int k, const string& seq,const unordered_set<minimizer>& A){
	minimizer kmer;
	double inter(0);

	for(unsigned int i(0);i+k<seq.size();++i){
		kmer=seq2int(seq.substr(i,k));
		if(A.count(kmer)>0){
			++inter;
		}
	}
	return double(100*inter/seq.size());
}




unordered_set<minimizer> filterUnitigs(const string& file, int k){
	ifstream in(file);
	unordered_set<minimizer> res;

	string line;
	while(!in.eof()){
		//~ getline(in,line);
		getline(in,line);
		//~ transform(line.begin(), line.end(),line.begin(), ::tolower);
		line=line.substr(0, line.size()-1);
		for(uint i(0);i+k<line.size();++i){
			res.insert(seq2int(line.substr(i,k)));
		}
	}
	return res;
}




void computeMinHash(size_t H, size_t k, size_t part, const vector<string>& reads, unordered_set<minimizer>* set, size_t L, size_t R){
	for (size_t i(L); i<R; ++i){
		string read(reads[i]);
		if(read.size()>100){
			vector<minimizer> sketch=minHashpart(H,k,read,part);
			myMutex.lock();
			for(size_t i(0);i<H;++i){
				set->insert(sketch[i]);
			}
			myMutex.unlock();
		}
	}
}




vector<size_t> bounds(size_t n,size_t size){
	vector<size_t> res;
	res.push_back(0);
	size_t d(size/n);

	for(size_t i(0); i<n-1;++i){
		res.push_back(i*d);
	}
	res.push_back(size);
	return res;
}




unordered_set<minimizer> filterUnitigs2(const vector<string>& V, size_t k, size_t H, size_t part){
	size_t nbThreads(8);
	unordered_set<minimizer> res;
	string line;
	vector<thread> threads;

	vector<size_t> limits = bounds(nbThreads, V.size());

	for (size_t i = 0; i < nbThreads; ++i) {
		threads.push_back(thread(computeMinHash,H,k,part,V,&res,limits[i],limits[i+1]));
	}

	for(auto &t : threads){t.join();}
	return res;
}




unordered_set<minimizer> filterContigsd(const string& file, int k){
	ifstream in(file);
	unordered_set<minimizer> res;

	string line;
	while(!in.eof()){
		getline(in,line);
		getline(in,line);
		transform(line.begin(), line.end(),line.begin(), ::tolower);
		for(uint i(0);i+k<line.size();++i){
			res.insert(seq2int(line.substr(i,k)));
		}
	}
	return res;
}




void indexSeqAux(const vector<string>& seqs, size_t H, size_t k, size_t part,  const unordered_set<minimizer>& filter, unordered_map<minimizer,unordered_set<uint32_t>>* index, uint32_t L, size_t R){

	for (uint32_t i(L); i<R; ++i){
		string seq=seqs[i];
		vector<minimizer> sketch;
		if(seq.size()<=(uint)H){
			sketch=allHash(k,seq);
		}else{
			sketch=minHashpart2(H,k,seq,part,filter);
			//~ sketch=minHash(H,k,unitig);
		}
		myMutex.lock();
		for(size_t j(0);j<sketch.size();++j){
			(*index)[sketch[j]].insert(i);
		}
		myMutex.unlock();
	}
}



unordered_map<minimizer,unordered_set<uint32_t>> indexSeq(const vector<string>& seqs, size_t H, size_t k, size_t part, const unordered_set<minimizer>& filter){
	size_t nbThreads(8);
	vector<thread> threads;
	unordered_map<minimizer,unordered_set<uint32_t>> index;
	vector<size_t> limits = bounds(nbThreads, seqs.size());

	for (size_t i(0); i<nbThreads; ++i) {
		threads.push_back(thread(indexSeqAux, seqs, H, k, part, filter, &index, limits[i], limits[i+1]));
	}

	for(auto &t : threads){t.join();}
	return index;
}




unordered_map<minimizer,unordered_set<uint32_t>> indexSeq2(const vector<string>& unitigs, int H, int k,int part){
	unordered_map<minimizer,unordered_set<uint32_t>> index;

#pragma omp parallel for
	for(uint i=0;i<unitigs.size();++i){
		string unitig=unitigs[i];
		vector<minimizer> sketch;
		if(unitig.size()<=(uint)H){
			sketch=allHash(k,unitig);
		}else{
			sketch=minHashpart(H,k,unitig,part);
			//~ sketch=minHash(H,k,unitig);
		}
		for(unsigned int j(0);j<sketch.size();++j){
#pragma omp critical
			{
				index[sketch[j]].insert(i);
			}
		}
	}
	return index;
}




vector<string> loadUnitigs(const string& unitigFile,bool homo){
	ifstream in(unitigFile);
	vector<string> res;
	int number(0);
	int size(0);

	string line;
	while(!in.eof()){
		getline(in,line);
		line=line.substr(0,line.size()-1);
		if(homo){
			line=homocompression(line);
		}
		if(line.size()>100){
			res.push_back(line);
			++number;
			size+=line.size();
		}
	}
	random_shuffle (res.begin(),res.end());
	cout<<"number of unitig : "<<number<<endl;
	cout<<"Mean size  : "<<size/number<<endl;
	return res;
}




unordered_map<string,vector<uint32_t>> getGraph(const vector<string>& unitigs,int k){
	unordered_map<string,vector<uint32_t>> graph;
	string unitig,seq1,seq2;

	for(uint i=0;i<unitigs.size();++i){
		unitig=unitigs[i];
		//~ cout<<unitig<<endl;
		seq1=getRepresent(unitig.substr(0,k));
		seq2=getRepresent(unitig.substr(unitig.size()-k,k));
		graph[seq1].push_back(i);
		if(seq1!=seq2){
			graph[seq2].push_back(i);
		}
	}

	return graph;
}




vector<string> loadContigs(const string& unitigFile){
	ifstream in(unitigFile);
	vector<string> res;

	string line;
	while(!in.eof()){
		getline(in,line);
		getline(in,line);
		//~ line=line.substr(0,line.size()-1);
		//~ string linerc(reversecompletment(line));
		transform(line.begin(), line.end(),line.begin(), ::tolower);
		//~ transform(linerc.begin(), linerc.end(),linerc.begin(), ::toupper);
		//~ cout<<line<<endl;
		//~ cin.get();
		res.push_back(line);
		//~ res.push_back(linerc);
		//~ res.push_back(reversecompletment(line.substr(0,line.size()-1)));
	}
	return res;
}




vector<string> loadFASTQ(const string& unitigFile,bool homo){
	ifstream in(unitigFile);
	vector<string> res;
	uint64_t n(0);
	uint64_t size(0);

	string line,lol;
	while(!in.eof()){		getline(in,lol);
		getline(in,line);
		getline(in,lol);
		getline(in,lol);
		transform(line.begin(), line.end(),line.begin(), ::tolower);
		if(line.size()>1000){
			if(homo){
				res.push_back(homocompression(line));
			}else{
				res.push_back(line);
			}
			++n;
			size+=line.size();
		}
	}
	random_shuffle (res.begin(),res.end());
	cout<<"number of reads : "<<n<<endl;
	cout<<"Mean  size : "<<size/n<<endl;
	return res;
}




void readFASTQ(const string& fastqFile, unordered_map<minimizer,unordered_set<uint32_t>>& index ,int H, int k, vector<string>& unitigs,int part,int multi){
	ifstream in(fastqFile);
	//~ int none(0),one(0),totalhit(0);
	int hit(0);
	int hitone(0);
	int reads(0);
	int goodread(0);
	int moyenread(0);
	int triplematchread(0);
	double coverage(0);
	cout<<"number of part:"<<part<<" number of matches required "<<multi<<endl;

#pragma omp parallel num_threads(8)
	{
		while(!in.eof()){
			double unitigLength(0);
			bool good=false;
			bool moyen=false;
			bool triplematch=false;
			string line,read;
			unordered_map<uint32_t,int> count;
#pragma omp critical
			{
				getline(in,line);
				getline(in,read);
				getline(in,line);
				getline(in,line);
			}
			if(read.size()>1000){
				transform(read.begin(), read.end(),read.begin(), ::tolower);
				//~ cout<<read<<endl;
				vector<minimizer> sketch;
				sketch=minHashpart(H,k,read,part);
				//~ sketch=minHash(H,k,read);
				//~ int hit(0);
				for(unsigned int j(0);j<sketch.size();++j){

					if(index.count(sketch[j])!=0){
						unordered_set<uint32_t> myset=index[(sketch[j])];
						for(auto it=myset.begin();it!=myset.end();++it){
							//~ cout<<unitigs[*it]<<endl;
							count[*it]++;
						}
					}

				}
#pragma omp atomic
				reads++;

				for(auto it=count.begin();it!=count.end();++it){
					if(it->second>=1){
#pragma omp atomic
						hitone++;
						moyen=true;
					}
					if(it->second>=multi){
						triplematch=true;
						//~ if(jaccard(12,unitigs[it->first],read)>0.05){
						if(true){
#pragma omp critical
							{
								hit++;
								unitigLength+=unitigs[it->first].size();
								//~ cout<<100*jaccard(k/2,unitigs[it->first],read)<<endl;
							}
							good=true;
						}
					}
				}
				if(good){
#pragma omp critical
					{
						goodread++;
						coverage+=(unitigLength/read.size());
						//~ cout<<(unitigLength/read.size())<<endl;
						//~ cout<<coverage<<endl;
					}
				}
				if(moyen){
#pragma omp atomic
					moyenread++;
				}
				if(triplematch){
#pragma omp atomic
					triplematchread++;
				}
				//~ cin.get();
			}
		}
	}
	cout<<"reads: "<<reads<<" 1hit: "<<hitone<<" 3hit: "<<hit<<" "<<100*hit/hitone<<" 1match/read(*10): "<<10*hitone/reads<<" 3match/reads(*10): "<<10*hit/reads<< " %Good read: "<<100*goodread/reads
	<<" %Read with a 3match: "<<100*triplematchread/reads<<" %Read with a 1match: "<<100*moyenread/reads<<" Coverage of good read: "<<coverage/goodread<<endl<<endl;
}




void readUnitigs2(const vector<string>& unitigs, unordered_map<minimizer,unordered_set<rNumber>>& index, size_t k, vector<string>& reads, size_t multi, size_t H, size_t part, size_t k2, double minJacc, const unordered_map<string,vector<uint32_t>>& graph){
	//~ int none(0),one(0),totalhit(0);
	//int k3(19);
	int hit(0),hitone(0),nreads(0),goodread(0),moyenread(0),triplematchread(0),jacc(0),failjacc(0),wrong(0),bigmistake(0),notmapped(0),mapped(0),missjacc(0),win(0),bigunitig(0);
	uint u(0);
	unordered_map <rNumber,int> read2coverage;
	unordered_map <rNumber,vector<string>> read2unitigs;
	unordered_set <rNumber> readCovered;
	unordered_set <rNumber> readHit;
	cout<<" number of matches required "<<multi<<endl;

#pragma omp parallel for num_threads(12)
	for(u=0;u<unitigs.size();++u){
		int localhit(0);
		bool good(false),moyen(false),triplematch(false);
		string line,unitig,read;
		unordered_map<rNumber,int> count;
		unordered_set<minimizer> min;
		unordered_map<rNumber,unordered_set<minimizer>> read2min;

		unitig=unitigs[u];
		if(unitig.size()>100){
			//~ cin.get();
			//~ cout<<"u:"<<u<<endl;
			unordered_set<minimizer> A;

			//~ transform(unitig.begin(), unitig.end(),unitig.begin(), ::tolower);
			if((int)unitig.size()<=H){
				for(unsigned int j(0);j+k<unitig.size();++j){
					uint32_t seq=seq2int(unitig.substr(j,k));
					if(min.count(seq)==0){
						if(index.count(seq)!=0){
							min.insert(seq);
							unordered_set<rNumber> myset=index[seq];
							for(auto it=myset.begin();it!=myset.end();++it){
								count[*it]++;
								read2min[*it].insert(seq);
							}
						}
					}
				}
			}else{
				vector<minimizer> sketch;
				sketch=minHashpart(H,k,unitig,part);

				//For each minimizer of the unitig
				for(unsigned int j(0);j<sketch.size();++j){
					minimizer seq=sketch[j];
					if(min.count(seq)==0){
						min.insert(seq);
						if(index.count(seq)!=0){
							unordered_set<rNumber> myset=index[seq];
							for(auto it=myset.begin();it!=myset.end();++it){
								++count[*it];
								read2min[*it].insert(seq);
							}
						}
					}
				}
			}

#pragma omp atomic
			nreads++;

			//For each read that share minimizers with the unitig
			for(auto it=count.begin();it!=count.end();++it){
				uint32_t readNumber=it->first;
				//~ cout<<"r:"<<it->first<<endl;
				bool goodreadb=false;
				//~ if(it->second>=(int)unitig.size()/multi){
				if(it->second>=multi){
#pragma omp critical
					{
						hitone++;
						readHit.insert(readNumber);
					}
					triplematch=true;
					unordered_set<minimizer> setMin(read2min[readNumber]);
					read=reads[readNumber];
					for(uint i(0);i+k<read.size() and !goodreadb;++i){
						//~ for(uint i(0);i<read.size();++i){
						minimizer seq(seq2int(read.substr(i,k)));
						if(setMin.count(seq)!=0){
							string region;

							if(i>1*unitig.size()){
								region=read.substr(i-1*unitig.size(),2*unitig.size());
							}else{
								region=read.substr(0,2*unitig.size());
							}


#pragma omp atomic
							jacc++;
							if(A.size()==0){
								vector<minimizer> sketchUnitig=allHash(k2,unitig);
								for(unsigned int i(0);i<sketchUnitig.size();++i){
									A.insert(sketchUnitig[i]);
								}
							}
							//~ if(region.size()>=read.size()){
							if(/* DISABLES CODE */ (false)){
								if(jaccard4(k2,region,A)>minJacc){
									moyen=true;
									//~ if(score(unitig,region)){
									localhit++;
#pragma omp critical
									{
										//~ hit++;
										readCovered.insert(readNumber);
										read2coverage[readNumber]+=read.size();
										read2unitigs[readNumber].push_back(unitig);
									}
									good=true;
									goodreadb=true;
								}else{
									i+=unitig.size();
#pragma omp critical
									{
										wrong++;
									}
								}
							}else{
								if(jaccard3(k2,region,A)>minJacc){
									moyen=true;
									//~ if(score(unitig,region)){
									localhit++;
#pragma omp critical
									{
										hit++;
										readCovered.insert(readNumber);
										read2coverage[readNumber]+=unitig.size();
										read2unitigs[readNumber].push_back(unitig);
									}
									good=true;
									goodreadb=true;

									//~ string seq=getRepresent(unitig.substr(unitig.size()-k3,k3));
									//~ vector<uint32_t> neigboor(graph.at(seq));
									//~ for(uint j(0);j<neigboor.size();++j){
									//~ string nextUnitig(unitigs[neigboor[j]]);
									//~ if(getRepresent(nextUnitig.substr(0,k3))==seq){
									//~ region=read.substr(i,2*nextUnitig.size());
									//~ if(jaccard3(k2,region,A)>minJacc){
									//~ #pragma omp atomic
									//~ win++;
									//~ }
									//~ }
									//~ }

								}else{
									i+=unitig.size();
#pragma omp critical
									{
										wrong++;
									}
								}
							}
							//~ }else{
							//~ #pragma omp atomic
							//~ failjacc++;
							//~ if(score(unitig,region) or score(reversecompletment(unitig),region)){
							//~ #pragma omp atomic
							//~ missjacc++;
							//~ }
							//~ }
						}
					}
					if(!goodreadb){
						//~ if(score(unitig,read) or score(reversecompletment(unitig),read)){
						//~ #pragma omp atomic
						//~ bigmistake++;
						//~ }else{
#pragma omp atomic
						notmapped++;
						//~ }
					}else
					{
#pragma omp atomic
						mapped++;
					}
				}

			}
			if(good){
#pragma omp critical
				{
					goodread++;
				}
			}
			if(moyen){
#pragma omp atomic
				moyenread++;
			}
			if(triplematch){
#pragma omp atomic
				triplematchread++;
			}
		}
	}
	cout<<"unitigs: "<<nreads<<" "<<unitigs.size()<<" hit: "<<hitone<<" aligns: "<<hit<<" wrong alignement: "<<wrong<<" %wrong alignment :"<<100*wrong/(hit+wrong)<<" match/unitigs: "<<hitone/(nreads)<<" align/unitigs: "<<hit/(nreads)<< " %unitig with an align: "<<100*goodread/nreads<<" %unitig with a match: "<<100*triplematchread/nreads<<" jacc:"<<jacc<<" jacc failed: "<<failjacc<<" %no jacc: "<<100*failjacc/(jacc)<<"%jacc fail: "<<100*missjacc/(failjacc+1)<<" %unitigs with a jacc "<<100*moyenread/nreads<<" read covered: "<<readCovered.size()<<" %read covered: "<<(100*readCovered.size())/reads.size()<<" %bigmistake: "<<(100*bigmistake)/(mapped+bigmistake)<<" "<<bigmistake<<" "<<mapped<<" "<<" %unmapped but hited reads "<<(100*notmapped)/(mapped+bigmistake+notmapped)<<" %read hit: "<<(100*readHit.size())/reads.size()<<" align/reads : "<<hit/reads.size()<<" win : "<<win<<" big unitig case : "<<bigunitig<<endl<<endl;
	uint64_t covered(0);
	uint64_t total(0);
	//~ for(auto it=read2unitigs.begin();it!=read2unitigs.end();++it){
	//~ cout<<"R :"<<reads[it->first]<<endl;
	//~ cout<<it->second.size()<<endl;
	//~ for(uint i(0);i<it->second.size();++i){
	//~ cout<<"U : "<<it->second[i]<<endl;
	//~ }
	//~ cout<<endl<<endl;
	//~ cin.get();
	//~ }

	for(auto it=read2coverage.begin();it!=read2coverage.end();++it){
		total+=reads[it->first].size();
		covered+=it->second;

	}
	cout<<"percent covered read: "<<covered/1000<<" "<<total/1000<<" "<<100*covered/total<<" "<<reads.size()<<endl;
	int l(0),n(0);
	for(uint i(0);i<reads.size();++i){
		if(readCovered.count(i)==0){
			l+=reads[i].size();
			++n;
		}
	}
	cout<<l<<" "<<n<<endl;
}


void readUnitigsAux(const vector<string>& unitigs, unordered_map<minimizer, unordered_set<rNumber>>* index, size_t k, const vector<string>& reads, size_t multi, size_t H, size_t part, size_t k2, double minJacc, const unordered_map<string,vector<uint32_t>>& graph, size_t L, size_t R){

	for (size_t i(L); i<R; ++i){
		string unitig=unitigs[i];
		unordered_set<minimizer> min;
		unordered_map<rNumber,size_t> count;
		unordered_map<rNumber,unordered_set<minimizer>> read2min;

		if(unitig.size()>100){
			unordered_set<minimizer> A;
			if(unitig.size()<=H){
				for(size_t j(0);j+k<unitig.size();++j){
					uint32_t seq=seq2int(unitig.substr(j,k));
					if(min.count(seq)==0){
						if(index->count(seq)!=0){
							min.insert(seq);
							unordered_set<rNumber> myset=index->at(seq);
							for(auto it=myset.begin();it!=myset.end();++it){
								count[*it]++;
								read2min[*it].insert(seq);
							}
						}
					}
				}
			}else{
				vector<minimizer> sketch;
				sketch=minHashpart(H,k,unitig,part);

				//For each minimizer of the unitig
				for(unsigned int j(0);j<sketch.size();++j){
					minimizer seq=sketch[j];
					if(min.count(seq)==0){
						min.insert(seq);
						if(index->count(seq)!=0){
							unordered_set<rNumber> myset=index->at(seq);
							for(auto it=myset.begin();it!=myset.end();++it){
								++count[*it];
								read2min[*it].insert(seq);
							}
						}
					}
				}
			}
			for(auto it=count.begin();it!=count.end();++it){
				uint32_t readNumber=it->first;
				bool goodreadb=false;
				if(it->second>=multi){
					unordered_set<minimizer> setMin(read2min[readNumber]);
					string read=reads[readNumber];
					for(uint i(0);i+k<read.size() and !goodreadb;++i){
						minimizer seq(seq2int(read.substr(i,k)));
						if(setMin.count(seq)!=0){
							string region;

							if(i>1*unitig.size()){
								region=read.substr(i-1*unitig.size(),2*unitig.size());
							}else{
								region=read.substr(0,2*unitig.size());
							}
							if(A.size()==0){
								vector<minimizer> sketchUnitig=allHash(k2,unitig);
								for(unsigned int i(0);i<sketchUnitig.size();++i){
									A.insert(sketchUnitig[i]);
								}
							}
							if(jaccard3(k2,region,A)>minJacc){
								//end
							}
						}
					}
				}
			}
		}
	}
	
}


void readUnitigs(const vector<string>& unitigs, unordered_map<minimizer,unordered_set<rNumber>>& index, size_t k, const vector<string>& reads, size_t multi, size_t H, size_t part, size_t k2, double minJacc, const unordered_map<string,vector<uint32_t>>& graph){
	size_t nbThreads(8);
	vector<thread> threads;
	vector<size_t> limits = bounds(nbThreads, unitigs.size());

	for (size_t i(0); i<nbThreads; ++i) {
		threads.push_back(thread(readUnitigsAux, unitigs, &index, k, reads, multi, H, part, k2, minJacc, graph, limits[i], limits[i+1]));
	}

	for(auto &t : threads){t.join();}
}






void readUnitigs2(const vector<string>& unitigs, unordered_map<minimizer, unordered_set<rNumber>>& index, size_t k, const vector<string>& reads, size_t multi, size_t H, size_t part, size_t k2, double minJacc, const unordered_map<string,vector<uint32_t>>& graph){
	size_t nbThreads(8);
	vector<thread> threads;
	vector<size_t> limits = bounds(nbThreads, unitigs.size());
	atomic<size_t> numberReads(0);

	for (size_t i(0); i<nbThreads; ++i) {
		threads.push_back(thread(readUnitigsAux,unitigs, &index, k, reads, multi, H, part, k2, minJacc, graph,limits[i], limits[i+1]));
	}

	for(auto &t : threads){t.join();}

}




void readContigs(const string& fastqFile, unordered_map<minimizer,unordered_set<uint32_t>>& index, int k, vector<string>& unitigs,int multi,int H,int part){
	ifstream in(fastqFile);
	//~ int none(0),one(0),totalhit(0);
	int hit(0);
	int hitone(0);
	int reads(0);
	int goodread(0);
	int moyenread(0);
	int triplematchread(0);
	double coverage(0);
	cout<<" number of matches required "<<multi<<endl;

#pragma omp parallel num_threads(8)
	{
		while(!in.eof()){
			double unitigLength(0);
			bool good=false;
			bool moyen=false;
			bool triplematch=false;
			string line,unitig;
			unordered_map<uint32_t,int> count;
			unordered_set<uint32_t> min;
#pragma omp critical
			{
				getline(in,unitig);
				getline(in,unitig);
			}
			//~ cout<<unitig<<endl;
			if(unitig.size()>3){
				transform(unitig.begin(), unitig.end(),unitig.begin(), ::tolower);

				if((int)unitig.size()<=H){

					for(unsigned int j(0);j+k<unitig.size();++j){
						uint32_t seq=seq2int(unitig.substr(j,k));
						if(min.count(seq)==0){
							if(index.count(seq)!=0){
								min.insert(seq);
								unordered_set<uint32_t> myset=index[seq];
								for(auto it=myset.begin();it!=myset.end();++it){
									count[*it]++;
								}
							}
						}
					}
				}else{

					vector<minimizer> sketch;
					sketch=minHashpart(H,k,unitig,part);

					for(unsigned int j(0);j<sketch.size();++j){
						uint32_t seq=sketch[j];
						if(min.count(seq)==0){
							if(index.count(seq)!=0){
								min.insert(seq);
								unordered_set<uint32_t> myset=index[seq];
								for(auto it=myset.begin();it!=myset.end();++it){
									count[*it]++;
								}
							}
						}
					}
				}

#pragma omp atomic
				reads++;

				for(auto it=count.begin();it!=count.end();++it){
					if(it->second>=1){
#pragma omp atomic
						hitone++;
						moyen=true;
					}
					if(it->second>=multi){
						triplematch=true;
						//~ if(jaccard(12,unitig,unitigs[it->first])>0.1){
						if(true){
#pragma omp critical
							{
								hit++;
								unitigLength+=unitigs[it->first].size();
								//~ cout<<100*jaccard(k/2,unitigs[it->first],read)<<endl;
							}
							good=true;
						}
					}
				}
				if(good){
#pragma omp critical
					{
						goodread++;
						coverage+=(unitigLength/unitig.size());
						//~ cout<<(unitigLength/read.size())<<endl;
						//~ cout<<coverage<<endl;
					}
				}
				if(moyen){
#pragma omp atomic
					moyenread++;
				}
				if(triplematch){
#pragma omp atomic
					triplematchread++;
				}
				//~ cout<<reads<<endl;
			}
		}
	}
	cout<<"unitigs: "<<reads<<" hit: "<<hitone<<" align: "<<hit<<" "<<100*hit/hitone<<" match/unitigs: "<<hitone/reads<<" align/unitigs: "<<hit/reads<< " %Good unitig: "<<100*goodread/reads
	<<" %unitig with a 3match: "<<100*triplematchread/reads<<" %unitig with a 1match: "<<100*moyenread/reads<<" Coverage of good unitig: "<<coverage/goodread<<endl<<endl;
}




void readFASTQforstat(const string& fastqFile,int H, int k){
	ifstream in(fastqFile);
	int max(0);
	int nread(0);
	string maxmin;
	unordered_map<minimizer,int> count;

#pragma omp parallel num_threads(8)
	{
		while(!in.eof()){
			string line,read;

#pragma omp critical
			{
				getline(in,line);
				getline(in,read);
				//~ getline(in,line);
				//~ getline(in,line);
				nread++;
			}
			if(read.size()>2){

				transform(read.begin(), read.end(),read.begin(), ::tolower);
				vector<minimizer> sketch=minHash(H,k,read);
				for(unsigned int j(0);j<sketch.size();++j){
#pragma omp critical
					{
						count[sketch[j]]++;
					}
				}
			}
		}
	}

	int diffmin(0);
	int nbmul(0);
	for(auto it=count.begin();it!=count.end();++it){
		diffmin++;
		if(it->second>2){
			nbmul++;
		}
		if(it->second>max){
			max=it->second;
			maxmin=it->first;
		}
	}
	cout<<max<<" "<<maxmin<<" "<<diffmin<<" "<<nbmul<<endl;
}




string compact(const string& seq1,const string& seq2,const string& kmer){
	size_t k(kmer.size()),s1(seq1.size()),s2(seq2.size());
	string rc2(reversecompletment(seq2));
	string rc1(reversecompletment(seq1));
	//~ cout<<"go "<<kmer<<endl;
	//~ cout<<seq1<<" "<<seq2<<endl;
	//~ cout<<reversecompletment(seq1)<<" "<<reversecompletment(seq2)<<endl;


	if(seq1.substr(0,k)==kmer){
		if(seq2.substr(s2-k,k)==kmer){
			return seq2+seq1.substr(k);
		}else{
			if(rc2.substr(s2-k,k)==kmer){
				return rc2+seq1.substr(k);
			}
			//~ else{return "";}
		}
	}
	if(seq2.substr(0,k)==kmer){
		if(seq1.substr(s1-k,k)==kmer){
			return seq1+seq2.substr(k);
		}else{
			if(rc1.substr(s1-k,k)==kmer){
				return rc1+seq2.substr(k);
			}
			//~ else{return "";}
		}
	}


	if(rc1.substr(0,k)==kmer){
		if(seq2.substr(s2-k,k)==kmer){
			return seq2+rc1.substr(k);
		}else{
			if(rc2.substr(s2-k,k)==kmer){
				return rc2+rc1.substr(k);
			}
			//~ else{return "";}
		}
	}

	if(rc2.substr(0,k)==kmer){
		if(seq1.substr(s1-k,k)==kmer){
			return seq1+rc2.substr(k);
		}else{
			if(rc1.substr(s1-k,k)==kmer){
				return rc1+rc2.substr(k);
			}
			//~ else{return "";}
		}
	}

	return "";
}




void readContigsforstats(const string& fastaFile,uint k, bool elag, bool compaction,bool unitigb){
	ifstream in(fastaFile);
	uint minSize(100);
	unordered_map<string,vector<int>> kmer2reads;
	int i(0);
	vector<string> unitigs;
	unordered_set<int> nottake;
	string line,read,seq1,seq2;

	while(!in.eof()){
		if(unitigb){
			//~ getline(in,line);
			getline(in,read);
			read=read.substr(0,read.size()-1);
		}else{
			getline(in,line);
			getline(in,read);
		}
		if(read.size()>2){
			transform(read.begin(), read.end(),read.begin(), ::tolower);
			unitigs.push_back(read);
			seq1=getRepresent(read.substr(0,k));
			seq2=getRepresent(read.substr(read.size()-k,k));
			kmer2reads[seq1].push_back(i);
			if(seq1!=seq2){
				kmer2reads[seq2].push_back(i);
			}
			i++;
		}
	}

	uint64_t length(0);
	int deadend(0),two(0),multiple(0),three(0),four(0),five(0),six(0),seven(0),eight(0),wtf(0),island(0);
	for(auto it=kmer2reads.begin();it!=kmer2reads.end();++it){

		if(it->second.size()==1){
			string unitig=unitigs[(it->second[0])];
			vector<int> v;
			if(unitig.size()>0){
				deadend++;
				if(elag){
					if(unitig.size()<=minSize){
						nottake.insert(it->second[0]);
					}
				}
				seq1=getRepresent(unitig.substr(0,k));
				seq2=getRepresent(unitig.substr(unitig.size()-k,k));

				if(seq1==it->first){
					v=kmer2reads[seq2];
				}else{
					if(seq2==it->first){
						v=kmer2reads[seq1];
					}else
					{
						wtf++;
					}
				}
				if(v.size()==1){
					island++;
					if(unitig.size()<=minSize and elag){
						nottake.insert(it->second[0]);
					}
				}
			}
		}

		if(it->second.size()==2){
			two++;
			if(compaction and (unitigs[it->second[0]].size()!=0) and (unitigs[(it->second[1])].size()!=0) ){
				string c(compact(unitigs[it->second[0]],unitigs[it->second[1]],it->first));
				if(c!=""){
					unitigs[it->second[0]]=c;
					unitigs[it->second[1]]="";
				}else{
					if((unitigs[it->second[0]]).size()<=minSize){
						unitigs[it->second[0]]="";
					}
					if((unitigs[it->second[1]]).size()<=minSize){
						unitigs[it->second[1]]="";
					}
				}
			}
		}

		if(it->second.size()==3){
			three++;
		}
		if(it->second.size()==4){
			four++;
		}
		if(it->second.size()==5){
			five++;
		}
		if(it->second.size()==6){
			six++;
		}
		if(it->second.size()==7){
			seven++;
		}
		if(it->second.size()==8){
			eight++;
		}
		if(it->second.size()>8){
			wtf++;
		}
		if(it->second.size()>1){
			multiple++;
		}
	}
	int n(0);
	//~ uint64_t length(0);
	ofstream out("unitig2.dot");
	for(uint i(0);i<unitigs.size();i++){
		if(nottake.count(i)==0 and unitigs[i].size()!=0){
			out<<unitigs[i]<<";"<<endl;
			length+=unitigs[i].size();
			n++;
		}
	}
	cout<<i<<" starting unitigs"<<endl;
	cout<<n<<" unitigs left"<<endl;
	cout<<length/n<<" mean length of unitigs"<<endl;
	cout<<i-n<<" unitigs removed"<<endl;
	cout<<island/2<<" island"<<endl;
	cout<<deadend<<" "<<two<<" "<<three<<" "<<four<<" "<<five<<" "<<six<<" "<<seven<<" "<<eight<<" "<<wtf<<" "<<multiple<<endl;
	cout<<endl;
}




unordered_set<string> getFilter(const string& faFile, int H,int k){
	ifstream in(faFile);
	unordered_set<string> minSet;

#pragma omp parallel num_threads(8)
	{
		while(!in.eof()){
			string line,read;
#pragma omp critical
			{
				getline(in,line);
				getline(in,read);
			}
			if(read.size()>3){
				//~ cout<<read<<endl;
				transform(read.begin(), read.end(),read.begin(), ::tolower);
				vector<string> sketch(minHashS(H,k,read));
				for(uint i(0);i<sketch.size();++i){
#pragma omp critical
					{
						minSet.insert(sketch[i]);
						//~ cout<<sketch[i]<<endl;
						//~ cin.get();
					}
				}
			}
		}
	}
	cout<<"size:"<<minSet.size()<<endl;
	return minSet;
}




unordered_set<string> getFilterTotal(const string& faFile,int k){
	ifstream in(faFile);
	unordered_set<string> minSet;

#pragma omp parallel num_threads(8)
	{
		while(!in.eof()){
			string line,read;
#pragma omp critical
			{
				getline(in,line);
				getline(in,read);
			}
			if(read.size()>3){
				transform(read.begin(), read.end(),read.begin(), ::tolower);
				vector<string> sketch(allHashS(k,read));
				for(uint i(0);i<sketch.size();++i){
#pragma omp critical
					{
						minSet.insert(sketch[i]);
					}
				}
			}
		}
	}
	return minSet;
}




int similar2filter(const string& faFile, int H,int k,int t, const unordered_set<string>& F){
	ifstream in(faFile);
	int res(0);
	unordered_set<string> minSet;

#pragma omp parallel num_threads(8)
	{
		while(!in.eof()){
			int n(0);
			string line,read;
#pragma omp critical
			{
				getline(in,line);
				getline(in,read);
			}
			if(read.size()>2){
				transform(read.begin(), read.end(),read.begin(), ::tolower);
				vector<string> sketch(minHashSf(H,k,read,F));
				for(uint i(0);i<sketch.size();++i){
					string min=sketch[i];
					if(F.count(min)!=0){
						++n;
					}
				}
				if(n>=t){
#pragma omp atomic
					++res;
				}
			}
		}
	}
	return res;
}




int similar2filterTotal(const string& faFile,int k,int t, const unordered_set<string>& B){
	ifstream in(faFile);
	int res(0);
	unordered_set<string> minSet;

#pragma omp parallel num_threads(8)
	{
		while(!in.eof()){
			int n(0);
			string line,read;
#pragma omp critical
			{
				getline(in,line);
				getline(in,read);
			}
			if(read.size()>2){
				transform(read.begin(), read.end(),read.begin(), ::tolower);
				vector<string> sketch(allHashS(k,read));
				for(uint i(0);i<sketch.size();++i){
					string min=sketch[i];
					if(B.count(min)!=0){
						++n;
						i+=k-1;
					}
				}
				if(n>=t){
#pragma omp atomic
					++res;
				}
			}
		}
	}
	return res;
}




int sizeFile(const string& faFile){
	ifstream in(faFile);
	int res(0);

#pragma omp parallel num_threads(8)
	{
		while(!in.eof()){
			string line,read;
#pragma omp critical
			{
				getline(in,line);
				getline(in,read);
			}
			if(read.size()>2){
#pragma omp atomic
				++res;
			}
		}
	}
	return res;
}




double similarity(const string& AFile,const string& BFile,int H, int k,int t){
	auto start=chrono::system_clock::now();
	auto Bfilter(getFilter(BFile,H,k));
	int AinB(similar2filter(AFile,H,k,t,Bfilter));
	cout<<AinB<<endl;
	auto Afilter(getFilter(AFile,H,k));
	int BinA(similar2filter(BFile,H,k,t,Afilter));
	cout<<BinA<<endl;
	//~ cout<<sizeFile(AFile)<<endl;
	//~ cout<<sizeFile(BFile)<<endl;
	double res((AinB+BinA)*100);
	res/=(sizeFile(AFile)+sizeFile(BFile));
	cout<<"similarity: "<<res<<endl;
	auto end1=chrono::system_clock::now();auto waitedFor=end1-start;
	cout<<"Done in "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;
	return res;
}




double similarityTotal(const string& AFile,const string& BFile, int k,int t){
	auto start=chrono::system_clock::now();
	auto Bfilter(getFilterTotal(BFile,k));
	int AinB(similar2filterTotal(AFile,k,t,Bfilter));
	cout<<AinB<<endl;
	auto Afilter(getFilterTotal(AFile,k));
	int BinA(similar2filterTotal(BFile,k,t,Afilter));
	cout<<BinA<<endl;
	//~ cout<<sizeFile(AFile)<<endl;
	//~ cout<<sizeFile(BFile)<<endl;
	double res((AinB+BinA)*100);
	res/=(sizeFile(AFile)+sizeFile(BFile));
	cout<<"similarity total: "<<res<<endl;
	auto end1=chrono::system_clock::now();auto waitedFor=end1-start;
	cout<<"Done in "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;
	return res;
}




int main(int argc, char ** argv) {
	size_t H1(100),k(15),part(1);
//	size_t k2(19);
	size_t k3(11);
	bool homo(false);
	srand((int)time(NULL));

	//~ similarityTotal("mini1.fasta","mini2.fasta",32,1);
	//~ similarity("mini1.fasta","mini2.fasta",10,32,1);
	//~
	//~ exit(0);

	//~ readContigsforstats("sim.contigs.fa",k2,false,false,false);

	//~ readContigsforstats("unitig.dot",k2,true,false,true);
	//~ for(int i(0);i<10;i++){
	//~ readContigsforstats("unitig2.dot",k2,true,false,true);
	//~ readContigsforstats("unitig2.dot",k2,false,true,true);
	//~ }

	//readContigsforstats("unitig2.dot",k2,false,false,true);

	auto start=chrono::system_clock::now();
	auto R(loadFASTQ("/Applications/PBMOG/Build/Products/Debug/PC10x_0001.fastq",homo));
	//~ auto R(loadFASTQ("PC20CE_0001.fastq",homo));
	auto U(loadUnitigs("/Applications/PBMOG/Build/Products/Debug/unitig2.dot",homo));
	//~ auto G(getGraph(U,k2));
	unordered_map<string,vector<uint32_t>> G;

	auto F(filterUnitigs2(U,k,H1,part));

	auto end1=chrono::system_clock::now();auto waitedFor=end1-start;
	cout<<"Reads loaded "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;

	auto index(indexSeq(R,H1,k,1*part,F));
	auto end2=chrono::system_clock::now();waitedFor=end2-end1;
	cout<<"Reads indexed "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;

	readUnitigs(U,index,k,R,2,H1,part,k3,100*(pow(1-2*0.10,k3)),G);
	auto end3=chrono::system_clock::now();waitedFor=end3-end2;
	cout<<"Unitigs mapped "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" seconds"<<endl<<endl;

	return 0;
}
