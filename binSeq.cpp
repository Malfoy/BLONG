//
//  binSeq.cpp
//  PBMOG
//
//  Created by malfoy on 27/01/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#include "binSeq.h"
#include <string>
#include <iostream>
#include <algorithm>



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
	0b10001101,0b10001001,0b10000101,0b10000001,0b10001100,0b10001000,0b10000100,0b10000000,
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




char revcompnadine (char s) {
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


//string revcomp (const string &s) {
//	string rc;
//	for (int i = (int)s.length() - 1; i >= 0; i--) rc += revcomp(s[i]);
//	return rc;
//}


string reversecomplementnadine (const string &s){
	string rc;
	for (int i = (int)s.length() - 1; i >= 0; i--) rc += revcompnadine(s[i]);
	return rc;

}






void printUC(unsigned char a){
	for (int i = 0; i < 8; i++) {
		printf("%d", ((a << i) & 0x80));
	}
	printf("\n");

}



unsigned char char2int( unsigned char c){
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



unsigned char int2char(unsigned char c){
	switch (c) {
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
	}
	if(c!=3){
		cout<<"wtf : ";
		printUC(c);
	}
	return 'T';
}




binSeq::binSeq(const string& str){
	isNumber=false;
	unsigned char mod(str.size()%3);
	unsigned char c;
	for (size_t i(0); i<str.size()-mod; i+=3){
		c=12;
		c+=char2int(str[i]);
		c<<=2;
		c+=char2int(str[i+1]);
		c<<=2;
		c+=char2int(str[i+2]);
		vect.push_back(c);
		//		printUC(c);
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
			c=1<<5;
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
	//	cout<<vect.size()<<endl;
	for(size_t i(0);i<vect.size();++i){
		unsigned char c(vect[i]);
		unsigned char mod(c/(1<<6));
		//		printUC(c);
		c<<=2;
		//		cout<<"mod : ";
		//		printUC(mod);

		switch (mod) {
			case 1:
				res.push_back(int2char(c/(1<<2)));
				break;
			case 2:
				res.push_back(int2char(c/(1<<4)));
				c<<=2;
				res.push_back(int2char((c/(1<<4))%4));
				break;

			case 3:
				res.push_back(int2char(c/(1<<6)));
				c<<=2;
				res.push_back(int2char(c/(1<<6)));
				c<<=2;
				res.push_back(int2char(c/(1<<6)));
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
		unsigned char c(vect[i]);
		unsigned char n(c/(1<<6));
		if(count+n<begin){
			count+=n;
		}else{
			if(count+n==begin){
				go=false;
			}else{
				go=false;
				unsigned char toGet((unsigned char)(count+n-begin));
				res.vect.push_back((unsigned char)(toGet<<6)+c%(1<<(6-2*(3-toGet))));
			}
		}
	}
	res.vect.insert(res.vect.end(), vect.begin()+i, vect.end());

	return res;
}


binSeq binSeq::sub(size_t begin, size_t size){
	binSeq res;
	size_t i(0);
	size_t countBegin(0);
	bool go(false);
	for (size_t count(0); count<size;) {
		unsigned char ch(vect[i]);
		unsigned char mod(ch/(1<<6));

		if(!go){
			if(countBegin+mod<begin){
				countBegin+=mod;
			}else{
				if(countBegin+mod==begin){
					go=true;
				}else{
					go=true;
					unsigned char toGet((unsigned char)(countBegin+mod-begin));
					res.vect.push_back((unsigned char)(toGet<<6)+ch%(1<<(6-2*(3-toGet))));
					count+=toGet;
				}
			}
		}else{
			if(count+mod<=size){
				res.vect.push_back(ch);
				count+=mod;
				//				cout<<1<<endl;
			}else{
				//				cout<<2<<endl;
				unsigned char toGet((unsigned char)(size-count));
				//				printUC(toGet);
				ch<<=2;
				res.vect.push_back((unsigned char)(toGet<<6)+ch/(1<<(8-2*(toGet))));
				//				printUC((unsigned char)(toGet<<6)+ch%(1<<(6-2*(3-toGet))));
				count+=toGet;
			}
		}

		++i;
	}

	return res;
}


/*
 uint64_t binSeq::getBegin(size_t size){
	uint64_t res(0);
	size_t i(0);
	for(size_t c(0); c<size;){
 unsigned char ch(vect[i]);
 unsigned char mod(ch/(1<<6));
 int64_t n(c+mod-size);
 if(n<=0){
 res<<=(2*mod);
 res+=ch%(1<<6);
 i++;
 c+=mod;
 }else{
 if(n==2){
 res<<=2;
 ch<<=2;
 res+=ch/(1<<6);
 c++;
 i++;
 }else{
 res<<=4;
 ch<<=2;
 res+=ch/(1<<4);
 c+=2;
 i++;
 }
 }
	}
	return res;
 }
 */


binSeq binSeq::getBegin(size_t size){
	binSeq res;
	size_t i(0);
	for(size_t c(0); c<size;){
		unsigned char ch(vect[i]);
		unsigned char mod(ch/(1<<6));
		int64_t n(c+mod-size);
		if(n<=0){
			//			res<<=(2*mod);
			//			res+=ch%(1<<6);
			res.vect.push_back(ch);
			i++;
			c+=mod;
		}else{
			if(n==2){
				//				res<<=2;
				ch<<=2;
				//				res+=ch/(1<<6);
				res.vect.push_back((1<<6)+ch/(1<<6));
				c++;
				i++;
			}else{
				//				res<<=4;
				ch<<=2;
				//				res+=ch/(1<<4);
				res.vect.push_back((2<<6)+ch/(1<<4));
				c+=2;
				i++;
			}
		}
	}
	return res;
}


/*

 uint64_t binSeq::getEnd(size_t size){
	uint64_t res(0);
	size_t i(vect.size()-1);
	for(size_t c(0); c<size;){
 unsigned char ch(vect[i]);
 unsigned char mod(ch/(1<<6));
 int64_t n(c+mod-size);
 if(n<=0){
 res+=(ch%(1<<6)<<(2*c));
 i--;
 c+=mod;
 }else{
 if(n==2){
 res+=(ch%(1<<2)<<(2*c));
 c++;
 i--;
 }else{
 res+=(ch%(1<<4)<<(2*c));
 c+=2;
 i--;
 }
 }
	}
	return res;
 }
 */

binSeq binSeq::getEnd(size_t size){
	binSeq res;
	size_t i(vect.size()-1);
	for(size_t c(0); c<size;){
		unsigned char ch(vect[i]);
		unsigned char mod(ch/(1<<6));
		int64_t n(c+mod-size);
		if(n<=0){
			//			res+=(ch%(1<<6)<<(2*c));
			res.vect.push_back(ch);
			i--;
			c+=mod;
		}else{
			if(n==2){
				//				res+=(ch%(1<<2)<<(2*c));
				res.vect.push_back((1<<6)+ch%(1<<2));
				c++;
				i--;
			}else{
				//				res+=(ch%(1<<4)<<(2*c));
				res.vect.push_back((2<<6)+ch%(1<<4));
				c+=2;
				i--;
			}
		}
	}
	::reverse(res.vect.begin(),res.vect.end());
	return res;
}



void binSeq::add(binSeq bs){
	vect.insert(vect.end(), bs.vect.begin(), bs.vect.end());
}


void binSeq::reverse(){
	vector<unsigned char> V;
	for(int i((int)vect.size()-1); i>-1; --i){
		V.push_back(rc[vect[i]]);
	}
	vect=V;
}


binSeq::binSeq(){

}


size_t binSeq::size(){
	return vect.size();
}


binSeq::binSeq(uint32_t n){
	isNumber=true;
	while(n!=0){
		vect.push_back((unsigned char)n%(1<<8));
		n>>=8;
	}
}


uint32_t binSeq::getNumber(){
	uint32_t res(0);
	for(int i((int)vect.size()-1); i>-1; --i){
		unsigned char c(vect[i]);
		res<<=8;
		res+=c;
	}
	return res;
}


uint64_t binSeq::getInt(){
	uint64_t res(0);
	for(size_t i(0); i<vect.size(); ++i){
		unsigned char c(vect[i]);
		unsigned char mod(c/(1<<6));
		res<<=(2*mod);
		res+=c%(1<<6);
	}
	return res;
}


void binSeq::clear(){
	vect.clear();
}


binSeq binSeq::getReverse(){
	binSeq res;
	for(int i((int)vect.size()-1); i>-1; --i){
		res.vect.push_back(rc[vect[i]]);
	}
	return res;
}





void testBinSeq(){
	cout<<"test start"<<endl;
	string str("TACCCATGCTAGCTGCATGACTGCTGACTGCATGTCGACTGATCGTACGTCGTAGCTGACTATATATAGCGCGCTATAAATATATACACACAGAGAGATTGTGTGTGTGTGTCGCGCGCGACACATATATAGTGTGTGCCAACACATTATATACG");
	binSeq bs(str);
	//	cout<<str<<endl;
	string str2(bs.str());
	//	cout<<str2/<<endl;
	if(str==bs.str()){
		cout<<"string constructor work and string export work"<<endl;
	}

	binSeq bs2(bs);
	if(str==bs2.str()){
		cout<<"binseq constructor work"<<endl;
	}

	bs.add(bs2);
	if(str+str==bs.str()){
		//		cout<<bs.str()<<endl;
		cout<<"add work"<<endl;
	}

	//	printUC((unsigned char)(bs.getBegin(4)));
	//	printUC((unsigned char)(bs.getEnd(4)));
	//	cout<<bs.getBegin(bs.getBegin(5));
	for(int i(0);i<(int)str.size();++i)
	{
		if(bs2.sub(i).str()==str.substr(i)){
			//			cout<<"sub work"<<endl;
		}else{
			cout<<bs2.sub(i).str()<<endl;
			cout<<"BUG SUB !!!!!!!"<<endl;
			cin.get();
		}

	}
	cout<<"sub work"<<endl;
	bs2.reverse();
	if(bs2.str()==reversecomplementnadine(str)){
		cout<<"reverse work"<<endl;
		cout<<reversecomplementnadine(str)<<endl;
	}else{
		cout<<bs2.str()<<endl;
		cout<<reversecomplementnadine(str)<<endl;
	}
	bs2.reverse();

	if(bs2.str()==str){
		cout<<"double reverse work"<<endl;
	}
	int a(6),b(120);
	if(bs2.sub(a,b).str()==str.substr(a,b)){
		cout<<"substr work with 2 arg"<<endl;
	}else{
		cout<<str.substr(a,b)<<endl;
		cout<<bs2.sub(a,b).str()<<endl;
	}

	uint n(2);
	binSeq bs3(n);
	if(n==bs3.getNumber()){
		cout<<"getInt and int constroctor work"<<endl;
	}else{
		cout<<bs3.getNumber()<<endl;
	}

	binSeq bs4(str);
	binSeq bs5(str);

	if((bs4+bs5).str()==str+str){
		cout<<"overload + work"<<endl;
	}
}








































