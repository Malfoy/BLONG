//
//  Utils.cpp
//  PBMOG
//
//  Created by malfoy on 09/01/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#include "Utils.h"

using namespace std;


size_t random(size_t max){
	default_random_engine generator;
	uniform_int_distribution<size_t> distribution(1,max);
	return(distribution(generator));
}

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
	for (int i = (int)s.length() - 1; i >= 0; i--){
		rc += revcomp(s[i]);
	}
	return rc;

}

bool isNuc(char c){
	//	if((65<c and c<116)){
	//		cout<<c<<endl;
	//	}
	return (c>64 and c<116);
}


string getRepresent (const string& str){
	string rc(reversecomplement(str));
	if(rc<str){
		return rc;
	}else{
		return str;
	}

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
	cout<<c<<endl;
	cin.get();
	return 0;
}

char int2nuc(char i){
	switch (i){
		case 0:
			return 'A';
			break;
		case 1:
			return 'C';
			break;
		case 2:
			return 'G';
			break;
		case 3:
			return 'T';
			break;
		default:
			return 'X';
			break;
	}
}

uint64_t xorshift64(uint64_t x) {
	x ^= x >> 12; // a
	x ^= x << 25; // b
	x ^= x >> 27; // c
	return x * UINT64_C(2685821657736338717);
}


vector<minimizer> allHash(size_t k,const string& seq){
	vector<minimizer> sketch;
	minimizer kmerS(seq2intStranded((seq.substr(0,k))));
	minimizer kmerRC(seq2intStranded((reversecomplement(seq.substr(0,k)))));
	minimizer kmer(min(kmerRC,kmerS));
	size_t i(0);
	do{
		sketch.push_back(kmer);
		if(i+k<seq.size()){
			updateMinimizer(kmerS, seq[i+k], k);
			updateMinimizerRC(kmerRC, seq[i+k], k);
			kmer=min(kmerRC,kmerS);
		}else{
			return sketch;
		}
		++i;
	}while(true);
	return sketch;
}


unordered_set <minimizer> allKmerSet(size_t k,const string& seq){
	unordered_set<minimizer> sketch;
	for(size_t i(0);i+k<=seq.size();++i){
		sketch.insert(seq2int(seq.substr(i,k)));
	}
	return sketch;
}


double scoreFromAlignment2(const string& seq){
	size_t begin(seq.size()),end(0);
	size_t errors(0),temp(0);
	for(size_t i(0);i<seq.size();++i){
		if(isNuc(seq[i])){
			if(begin==seq.size()){
				begin=i;
			}
			errors+=temp;
			temp=0;
			end=i;
		}else{
			if(begin!=seq.size()){
				temp++;
			}
		}
	}
	//	cout<<errors<<" "<<begin<<" "<<end<<endl;
	double res((100*errors)/(end-begin));
	return res;
}

double scoreFromAlignment(const string& seq1,const string& seq2){
	size_t errors(0);
	for(size_t i(0);i<seq1.size();++i){
		if(seq1[i]!=seq2[i]){
			++errors;
		}
	}
	double res((100*errors)/(seq1.size()));
	return res;
}

unordered_set <minimizer> allKmerSetStranded(size_t k,const string& seq){
	unordered_set<minimizer> sketch;
	minimizer min(seq2intStranded(seq.substr(0,k)));
	size_t i(0);
	do{
		sketch.insert(min);
		if(i+k<seq.size()){
			updateMinimizer(min, seq[i+k], k);
		}else{
			return sketch;
		}
		++i;
	}while(true);
	return sketch;
}


unordered_multimap<string,string> allKmerMapStranded(size_t k,const string& seq, char nuc){
	unordered_multimap<string,string> sketch;
	for(size_t i(0); i+k<=seq.size(); ++i){
		string kmer(seq.substr(i,k));
		//		if(kmer.size()!=k){
		//			cout<<"watfw"<<endl;
		//		}
		sketch.insert({kmer.substr(0,nuc),kmer.substr(nuc)});
		//				cout<<kmer.substr(0,nuc)<<" "<<kmer.substr(nuc)<<endl;
	}

	return sketch;
}


vector<string> kmerCounting(const string& fileName,size_t k){
	ifstream in(fileName);
	unordered_set<string> set;
	vector<string> res;
	string read,line;
	getline(in,line);
	getline(in,line);
	ofstream out("kmers.dot");
	//	cout<<read<<endl;
	read+=line;
	while (in.peek()=='A' or in.peek()=='C' or in.peek()=='G' or in.peek()=='T') {
		getline(in,line);
		read+=line;
	}
	for(size_t i(0);i<read.size() ;++i){
		string kmer(getRepresent((read.substr(i,k))));
		if(kmer.size()==k){
			//		cout<<kmer<<endl;
			set.insert(kmer);
		}else{
			break;
		}
	}

	for (auto it=set.begin(); it!=set.end(); ++it){
		string str(*it);
		transform(str.begin(), str.end(), str.begin(), ::tolower);
		out<<str<<":"<<endl;
		//		res.push_back(*it);
	}


	return res;
}


double jaccard(size_t k, const string& seq,const unordered_set<minimizer>& genomicKmers){
	minimizer kmer;
	double inter(0);
	//	cout<<"jacc "<<seq<<k<<endl;

	for(size_t i(0);i+k<=seq.size();++i){
		kmer=seq2int(seq.substr(i,k));
		//		cout<<seq.substr(i,k)<<endl;;
		if(genomicKmers.unordered_set::count(kmer)>0){
			//			cout<<"hit"<<endl;
			++inter;
		}
	}
	//	return double(100*inter/(genomicKmers.size()));
	//	return double(100*inter/(seq.size()-k));
	//	cout<<max(double(100*inter/(genomicKmers.size())),double(100*inter/(seq.size()-k)))<<endl;
	return max(double(100*inter/(genomicKmers.size())),double(100*inter/(seq.size()-k)));
}


double jaccardStranded(size_t k, const string& seq, const unordered_set<minimizer>& genomicKmers){
	minimizer kmer(seq2intStranded(seq.substr(0,k)));
	double inter(0);
	size_t i(0);
	do{
		if(genomicKmers.unordered_set::count(kmer)>0){
			++inter;
		}
		if(i+k<seq.size()){
			updateMinimizer(kmer, seq[i+k], k);
		}else{
			//			cout<<inter<<" "<<genomicKmers.size()<<" "<<seq.size()<<endl;
			return double(100*inter/(seq.size()-k+1));
			//			return max(double(100*inter/(genomicKmers.size())),double(100*inter/(seq.size()-k)));
		}
		++i;
	}while(true);
	//	return double(100*inter/(genomicKmers.size()));
	//	return double(100*inter/(seq.size()-k));
	//	return max(double(100*inter/(genomicKmers.size())),double(100*inter/(seq.size()-k)));
}


bool equalStr(const string& seq1, const string& seq2){
	size_t size(min(seq1.size(),seq2.size()));
	return (seq1.substr(0,size))==seq2.substr(0,size);
}


bool isCorrect(const string& seq,const string& ref){
	for(size_t i(0); i<seq.size(); ++i){
		if(seq[i]!=ref[i]){
			if(seq[i+1]==ref[i]){
				return equalStr(seq.substr(i+2),ref.substr(i+1));
			}
			if(seq[i]==ref[i+1]){
				return equalStr(seq.substr(i+1),ref.substr(i+2));
			}
			return (seq.substr(i+1)==ref.substr(i+1));
		}
	}
	return true;
}


double jaccardStrandedErrors(size_t k, const string& seq, const unordered_multimap<string, string>& genomicKmers, char nuc){
	double inter(0);
	string kmer;
	kmer.reserve(k);
	size_t i(0);
	for(; i+k<=seq.size(); ++i){
		kmer=seq.substr(i,k);
		if(kmer.size()!=k){
			cout<<"wtf"<<endl;
		}
		auto range(genomicKmers.equal_range(kmer.substr(0,nuc)));
		for (auto it(range.first); it!=range.second; it++){
			if(isCorrect(kmer.substr(nuc),it->second)){
				//				if(isCorrect("ATCGTTTT", "ATTGTTTT")){
				//					cout<<"true"<<endl;
				//					cin.get();
				//				}else{
				//					cout<<"false"<<endl;
				//					cin.get();
				//				}
				//				if(kmer.substr(5)!=it->second){
				//					cout<<kmer.substr(5)<<" "<<it->second<<endl;
				//					cin.get();
				//				}
				inter++;
				break;
			}else{
				//				cout<<"fail : "<<kmer.substr(5)<<" "<<it->second<<endl;
			}
		}
	}
	//	cout<<i<<" "<<seq.size()-k+1<<endl;
	return double(100*inter/(seq.size()-k+1));;
}





void minHash2(size_t H, size_t k, const string& seq, vector<minimizer>& previous){
	vector<uint64_t> sketch(H);
	vector<minimizer> sketchs(H);
	uint64_t hashValue;
	//	hash<uint32_t> hash;

	minimizer kmerS(seq2intStranded(seq.substr(0,k)));
	minimizer kmerRC(seq2intStranded(reversecomplement(seq.substr(0,k))));
	minimizer kmer(min(kmerS,kmerRC));
	//	hashValue=hash(kmer);
	hashValue=xorshift64(kmer);
	for(size_t j(0); j<H; ++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=xorshift64(hashValue);
	}
	for(size_t i(1); i+k<seq.size(); ++i){
		updateMinimizer(kmerS, seq[i+k], k);
		updateMinimizerRC(kmerRC, seq[i+k], k);
		kmer=(min(kmerS,kmerRC));
		hashValue=xorshift64(kmer);
		//		hashValue=hash(kmer);
		for(size_t j(0); j<H; ++j){
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
	//	cout<<"lol"<<endl;
	//	cin.get();
	minimizer res(0);
	for(uint i(0);i<seq.size();++i){
		res<<=2;
		res+=nuc2int(str[i]);
	}
	return res;
}


minimizer seq2intStranded(const string& seq){
	minimizer res(0);
	for(uint i(0);i<seq.size();++i){
		res<<=2;
		res+=nuc2int(seq[i]);
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


string compaction(const string& seq1,const string& seq2, size_t k){
	size_t s1(seq1.size()),s2(seq2.size());
	if(s1==0){return "";}
	if(s2==0){return "";}

	string rc2(reversecomplement(seq2));
	//	string rc1(reversecomplement(seq1));

	string end1(seq1.substr(s1-k,k));
	string beg2(seq2.substr(0,k));
	if(end1==beg2){
		return seq1+(seq2.substr(k));
	}

	string begrc2(rc2.substr(0,k));
	if(end1==begrc2){
		return seq1+(rc2.substr(k));
	}

	string beg1(seq1.substr(0,k));
	string end2(seq2.substr(s2-k,k));
	if(beg1==end2){
		return seq2+(seq1.substr(k));
	}

	string endrc2(rc2.substr(s2-k,k));
	if(beg1==endrc2){
		return rc2+(seq1.substr(k));
	}
	//	cout<<"failt"<<endl;
	return "";
}

string compactionEnd(const string& seq1,const string& seq2, size_t k){

	size_t s1(seq1.size()),s2(seq2.size());
	if(s1==0 or s2==0){return "";}

	string rc2(reversecomplement(seq2));
	//		string rc1(reversecomplement(seq1));

	string end1(seq1.substr(s1-k,k));
	string beg2(seq2.substr(0,k));
	if(end1==beg2){
		return seq1+(seq2.substr(k));
	}

	string begrc2(rc2.substr(0,k));
	if(end1==begrc2){
		return seq1+(rc2.substr(k));
	}
	return "";
}

string compactionBegin(const string& seq1,const string& seq2, size_t k){
	size_t s1(seq1.size()),s2(seq2.size());
	if(s1==0 or s2==0){return "";}

	string rc2(reversecomplement(seq2));
	//	string rc1(reversecomplement(seq1));

	string beg1(seq1.substr(0,k));
	string end2(seq2.substr(s2-k,k));
	if(beg1==end2){
		return seq2+(seq1.substr(k));
	}

	string endrc2(rc2.substr(s2-k,k));
	if(beg1==endrc2){
		return rc2+(seq1.substr(k));
	}
	//	cout<<"failt"<<endl;
	return "";
}



void readContigsforstats(const string& File, size_t k, bool elag, bool compact,bool unitigb){
	ifstream in(File);
	uint minSize(50);
	unordered_map<string,vector<size_t>> kmer2reads;
//	kmer2reads.set_empty_key("0");
	int i(0);
	vector<string> unitigs;
	unordered_set<size_t> nottake;
	string line,read,seq1,seq2;

	while(!in.eof()){
		if(unitigb){
			//~ getline(in,line);
			getline(in,read);
			read=read.substr(0,read.size()-1);
			transform(read.begin(), read.end(),read.begin(), ::toupper);
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
//					if(unitig.size()<=minSize and elag){
//						nottake.insert(it->second[0]);
//					}
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
		if(nottake.unordered_set::count(ii)==0 and unitigs[ii].size()!=0){
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
//	graph.set_empty_key("0");
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


int positionInSeq(const string& seq, minimizer min, size_t k){
	minimizer kmer(seq2intStranded(seq.substr(0,k)));
	minimizer kmerRC(seq2intStranded(reversecomplement((seq.substr(0,k)))));
	for(int i(0);;++i){
		if(min==kmer or min==kmerRC){
			return i;
		}

		if(i+k<seq.size()){
			updateMinimizer(kmer, seq[i+k], k);
			updateMinimizerRC(kmerRC, seq[i+k], k);
		}else{
			return -1;
		}
	}
	return -1;
}


int positionInSeqStranded(const string& seq, minimizer min, size_t k){
	minimizer kmer(seq2intStranded(seq.substr(0,k)));
	for(int i(0);;++i){
		if(min==kmer){
			return (int)i;
		}

		if(i+k<seq.size()){
			updateMinimizer(kmer, seq[i+k], k);
		}else{

			return -1;
		}
	}
	return -1;
}


int positionInSeqStrandedEnd(const string& seq, minimizer min, size_t k){
	//	cout<<seq<<endl;
	minimizer kmer(seq2intStranded(seq.substr(seq.size()-k,k)));
	for(int i((int)seq.size()-(int)k);; --i){
		//		printMinimizer(kmer, k);
		if(min==kmer){
			return (int)i;
		}
		if(i>=0){
			updateMinimizerEnd(kmer, seq[i], k);
		}else{
			//			cin.get();
			return -1;
		}
	}
	return -1;
}





vector<string> loadFASTQ(const string& unitigFile,bool homo,size_t sizeMin,char frac){
	ifstream in(unitigFile);
	vector<string> res;
	int n(0),rn(0);
	uint64_t size(0);
	//TODO compatibility fastq and fasta plz
	string line,lol;
	while(!in.eof()){
		getline(in,lol);
		getline(in,line);
		if(in.peek()=='+'){
			getline(in,lol);
			getline(in,lol);
		}
		if(line.size()>sizeMin){
			if(homo){
				res.push_back(homocompression(line));
			}else{
				//				res.push_back(randomString(10000));
				if(n%frac==0){
					res.push_back(line);
					rn++;
				}
			}
			++n;
			size+=line.size();
		}
	}
	random_shuffle (res.begin(),res.end());
	cout<<"number of reads : "<<rn<<endl;
	cout<<"Mean  size : "<<size/(max(1,n))<<endl;
	return res;
}

vector<string> loadFASTA(const string& unitigFile,bool homo,size_t sizeMin, size_t frac){
	ifstream in(unitigFile);
	vector<string> res;
	int n(0);
	uint64_t size(0);
	//TODO compatibility fastq and fasta plz
	string line,lol,more;
	while(!in.eof()){
		getline(in,lol);
		getline(in,line);
		while(in.peek()!='>' and !in.eof()){
			getline(in,more);
			line+=more;
		}
		if(line.size()>sizeMin){
			++n;
			if(homo){
				res.push_back(homocompression(line));
			}else{
				//				res.push_back(randomString(10000));
				if(n%frac==0){
					res.push_back(line);
				}
			}
			size+=line.size();
		}

	}
	random_shuffle (res.begin(),res.end());
	cout<<"number of reads : "<<n<<endl;
	cout<<"Mean  size : "<<size/(max(1,n))<<endl;
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

string randomString( size_t length )
{
	auto randchar = []() -> char
	{
		const char charset[] ="ATCG";
		const size_t max_index = (sizeof(charset) - 1);
		return charset[ rand() % max_index ];
	};
	string str(length,0);
	generate_n( str.begin(), length, randchar );
	return str;
}

void printMinimizer(minimizer min,size_t k){
	string res;
	for(size_t i(0); i<k; ++i){
		res+=int2nuc(min%4);
		min>>=2;
	}
	reverse(res.begin(),res.end());
	cout<<res<<endl;
}

void updateMinimizer(minimizer&	min, char nuc,size_t k){
	minimizer offset(1<<(2*k));
	min<<=2;
	min+=nuc2int(nuc);
	min%=offset;
}

void updateMinimizerEnd(minimizer&	min, char nuc,size_t k){
	min>>=2;
	min+=(nuc2int(nuc)<<(2*k-2));
}

void updateMinimizerRC(minimizer&	min, char nuc,size_t k){
	//	printMinimizer(min, k);
	//	cout<<nuc<<endl;
	min>>=2;
	min+=((3-nuc2int(nuc))<<(2*k-2));
	//	printMinimizer(min, k);
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
			//			read=read.substr(2000+random(read.size()-4000),200);
			res.push_back(read);
			//			res.push_back("ACATCAAAGCTAGTGTGAGCTCCGATAATCACTGTGAGAAAAGGCGATAGGAACCGCATGACTCCAATGTAGGTCCTTCCCGGGTGGGGACCTGGCGTGAGGCAGACTGCGGCCGATGGTGAGAAAGGAATTCAATGAGTTGCATCGGCACCCCGAAGTATAGCACGTAGGTCAGGACGTTTCTTATACAGAGCCCGGAA");
			//			cout<<read<<endl;
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
	uint64_t hashValue;;
	//~ hash<uint32_t> hash;

	minimizer kmerS=seq2intStranded(seq.substr(0,k));
	minimizer kmerRC=seq2intStranded(reversecomplement(seq.substr(0,k)));
	minimizer kmer(min(kmerS,kmerRC));
	hashValue=xorshift64(kmer);
	//~ hashValue=hash(kmer);
	for(size_t j(0);j<H;++j){
		sketch[j]=hashValue;
		sketchs[j]=kmer;
		hashValue=xorshift64(hashValue);
	}

	for(size_t i(1);i+k<seq.size();++i){
		updateMinimizerRC(kmerRC, seq[i+k], k);
		updateMinimizer(kmerS, seq[i+k], k);
		kmer=min(kmerRC,kmerS);

		if(filter.unordered_set::count(kmer)!=0){
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



