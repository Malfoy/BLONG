//
//  Utils.cpp
//  PBMOG
//
//  Created by malfoy on 09/01/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#include "Utils.h"

using namespace std;


string reversecompletment2(const string& str){
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


char revcomp (char s) {
	if (s == 'A') return 'T';
	else if (s == 'C') return 'G';
	else if (s == 'G') return 'C';
	else if (s == 'T') return 'A';
	else if (s == 'a') return 't';
	else if (s == 'c') return 'g';
	else if (s == 'g') return 'c';
	else if (s == 't') return 'a';
	return 'X';
}


string revcomp (const string &s) {
	string rc;
	for (int i = (int)s.length() - 1; i >= 0; i--) rc += revcomp(s[i]);
	return rc;
}


string getRepresent2(const string &s){
	string rc;
	for (int i = 0; i < (int)s.length(); i++) {
		char c = revcomp(s[s.length() - 1 - i]);
		if (s[i] < c) {
			return s;
		} else if (s[i] > c) {
			return revcomp(s);
		}
	}
	return revcomp(s);
}


string reversecomplement (const string &s){
	string rc;
	for (int i = (int)s.length() - 1; i >= 0; i--) rc += revcomp(s[i]);
	return rc;

}


string getRepresent (const string& str){
	return(min(str,reversecomplement(str)));
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


unordered_set <minimizer> allKmerSet(size_t k,const string& seq){
	unordered_set<minimizer> sketch;
	for(size_t i(0);i+k<seq.size();++i){
		sketch.insert(seq2int(seq.substr(i,k)));
	}
	return sketch;
}


double jaccard(size_t k, const string& seq,const unordered_set<minimizer>& genomicKmers){
	minimizer kmer;
	double inter(0);

	for(size_t i(0);i+k<seq.size();++i){
		kmer=seq2int(seq.substr(i,k));
		if(genomicKmers.count(kmer)>0){
			++inter;
		}
	}
//	return double(100*inter/(genomicKmers.size()));
//	return double(100*inter/(seq.size()-k));
	return max(double(100*inter/(genomicKmers.size())),double(100*inter/(seq.size()-k)));
}


double jaccardAlt(size_t k, const string& seq,const unordered_set<minimizer>& genomicKmers){
	minimizer kmer;
	double inter(0);

	for(size_t i(0);i+k<seq.size();++i){
		kmer=seq2int(seq.substr(i,k));
		if(genomicKmers.count(kmer)>0){
			++inter;
		}
	}
	return double(100*inter/(seq.size()-k));
}


void minHash2(size_t H, size_t k, const string& seq, vector<minimizer>& previous){
	vector<uint64_t> sketch(H);
	vector<minimizer> sketchs(H);
	uint64_t hashValue;
	minimizer kmer;
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


minimizer seq2int(const string& seq){
	string str(getRepresent(seq));
	minimizer res(0);
	for(uint i(0);i<seq.size();++i){
		res+=nuc2int(str[i]);
		res<<=2;
	}

	return res;
}


vector<minimizer> minHashpart(size_t H, size_t k,const string& seq, size_t part){
	vector<minimizer> result;
	size_t size(seq.size()/part);
	for(size_t i(0);i<part;++i){
		minHash2(H/part,k,seq.substr(i*size,size+k),result);
	}
	return result;
}


vector<size_t> bounds(size_t n,size_t size){
	vector<size_t> res;
	res.push_back(0);
	size_t d(size/n);

	for(size_t i(1); i<n;++i){
		res.push_back(i*d);
	}
	res.push_back(size);
	return res;
}


string compaction(string& seq1, string& seq2, size_t k){
	size_t s1(seq1.size()),s2(seq2.size());
	if(s1==0){return seq2;}
	if(s2==0){return seq1;}

	string rc2(reversecomplement(seq2));
	string rc1(reversecomplement(seq1));


	if(seq1.substr(0,k)==seq2.substr(s2-k,k)){
		return seq2+seq1.substr(k);
	}else{
		if(rc2.substr(s2-k,k)==seq1.substr(0,k)){
			return rc2+seq1.substr(k);
		}
	}

	if(seq2.substr(0,k)==seq1.substr(s1-k,k)){
		return seq1+seq2.substr(k);
	}else{
		if(rc1.substr(s1-k,k)==seq2.substr(0,k)){
			return rc1+seq2.substr(k);
		}
	}

	if(rc1.substr(0,k)==seq2.substr(s2-k,k)){
		return seq2+rc1.substr(k);
	}else{
		if(rc2.substr(s2-k,k)==rc1.substr(0,k)){
			return rc2+rc1.substr(k);
		}
	}

	if(rc2.substr(0,k)==seq1.substr(s1-k,k)){
		return seq1+rc2.substr(k);
	}else{
		if(rc1.substr(s1-k,k)==rc2.substr(0,k)){
			return rc1+rc2.substr(k);
		}
	}
	return "";
}



void readContigsforstats(const string& File, size_t k, bool elag, bool compact,bool unitigb){
	ifstream in(File);
	uint minSize(100);
	unordered_map<string,vector<size_t>> kmer2reads;
	kmer2reads.set_empty_key("0");
	int i(0);
	vector<string> unitigs;
	unordered_set<size_t> nottake;
	string line,read,seq1,seq2;

	while(!in.eof()){
		if(unitigb){
			//~ getline(in,line);
			getline(in,read);
			read=read.substr(0,read.size()-1);
		}else{
			read="";
			getline(in,line);
			getline(in,line);
			read+=line;
			while (in.peek()=='A' or in.peek()=='C' or in.peek()=='G' or in.peek()=='T') {
				getline(in,line);
				read+=line;
			}
		}
		if(read.size()>2){
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
			vector<size_t> v;
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
			if(compact and (unitigs[it->second[0]].size()!=0) and (unitigs[(it->second[1])].size()!=0) ){
				string c(compaction( unitigs[it->second[0]], unitigs[it->second[1]], k));
				if(c!=""){
					unitigs[it->second[0]]=c;
					unitigs[it->second[1]]="";
				}else{
//					if((unitigs[it->second[0]]).size()<=minSize){
//						unitigs[it->second[0]]="";
//					}
//					if((unitigs[it->second[1]]).size()<=minSize){
//						unitigs[it->second[1]]="";
//					}
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
			cout<<"wut"<<endl;
		}
		if(it->second.size()>1){
			multiple++;
		}
	}
	int n(0);
	//~ uint64_t length(0);
	ofstream out("unitigClean.fa");
	for(uint ii(0);ii<unitigs.size();ii++){
		if(nottake.count(ii)==0 and unitigs[ii].size()!=0){
			out<<">"<<ii<<endl;
			out<<unitigs[ii]<<endl;
			length+=unitigs[ii].size();
			n++;
		}
	}
	cout<<i<<" starting unitigs"<<endl;
	cout<<n<<" unitigs left"<<endl;
	cout<<length/(n+1)<<" mean length of unitigs"<<endl;
	cout<<i-n<<" unitigs removed"<<endl;
	cout<<island/2<<" island"<<endl;
	cout<<deadend<<" "<<two<<" "<<three<<" "<<four<<" "<<five<<" "<<six<<" "<<seven<<" "<<eight<<" "<<wtf<<" "<<multiple<<endl;
	cout<<endl;
}



unordered_map<string,vector<uNumber>> getGraph(const vector<string>& unitigs, size_t k){
	unordered_map<string,vector<uNumber>> graph;
	graph.set_empty_key("0");
	string unitig,seq1,seq2;

	for(uint32_t i(0);i<unitigs.size();++i){
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



vector<string> loadFASTQ(const string& unitigFile,bool homo){
	ifstream in(unitigFile);
	vector<string> res;
	uint64_t n(0);
	uint64_t size(0);

	string line,lol;
	while(!in.eof()){
		getline(in,lol);
		getline(in,line);
		getline(in,lol);
		getline(in,lol);
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
	cout<<"Mean  size : "<<size/(n+1)<<endl;
	return res;
}


string homocompression(const string& seq){
	string res;
	res.push_back(seq[0]);
	char last(seq[0]);
	for(uint i(1);i<seq.size();++i){
		if(seq[i]!=last){
			res.push_back(seq[i]);
			last=seq[i];
		}
	}
	return res;
}


vector<string> loadUnitigs(const string& unitigFile,bool homo){
	ifstream in(unitigFile);
	vector<string> res;
	int number(0);
	int size(0);

	string line,read;
	while(!in.eof()){
		read="";
		getline(in,line);
		getline(in,line);
		read+=line;
		while (in.peek()=='A' or in.peek()=='C' or in.peek()=='G' or in.peek()=='T') {
			getline(in,line);
			read+=line;
		}
		if(homo){
			line=homocompression(read);
		}
		if(read.size()>2){
			res.push_back(read);
			++number;
			size+=read.size();
		}
	}
	random_shuffle (res.begin(),res.end());
	cout<<"number of unitig : "<<number<<endl;
	cout<<"Mean size of unitigs : "<<size/(number+1)<<endl;
	return res;
}




void minHash3(size_t H, size_t k,const string& seq, vector<minimizer>& previous, const unordered_set<minimizer>& filter){
	vector<uint64_t> sketch(H);
	vector<minimizer> sketchs(H);
	uint64_t hashValue;
	minimizer kmer;
	//~ hash<uint32_t> hash;

	kmer=seq2int(seq.substr(0,k));
	hashValue=xorshift64(kmer);
	//~ hashValue=hash(kmer);
	for(size_t j(0);j<H;++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=xorshift64(hashValue);
	}

	for(size_t i(1);i+k<seq.size();++i){
		kmer=seq2int(seq.substr(i,k));
		if(filter.count(kmer)!=0){
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
	}
	previous.insert(previous.end(),sketchs.begin(),sketchs.end());
}


uint64_t xorshift(uint64_t x) {
	x ^= x >> 12; // a
	x ^= x << 25; // b
	x ^= x >> 27; // c
	return x * UINT64_C(2685821657736338717);
}








vector<minimizer> minHashpart2(size_t H, size_t k,const string& seq, size_t part, const unordered_set<minimizer>& filter){
	vector<minimizer> result;
	size_t size(seq.size()/part);
	for(size_t i(0);i<part;++i){
		//~ minHash3(H/part,k,seq.substr(i*size,size+16),result,filter);
		minHash3(H/part,k,seq.substr(i*size,size),result,filter);
	}
	return result;
}



