#include<stdio.h>
#include<random>
#include<iostream>
#include<fstream>
using namespace std;

void password(){

	ofstream f("password.dat");

	random_device rnd;
	mt19937 mt(rnd());
	uniform_int_distribution<> prior_dist(0,25);
	uniform_int_distribution<> pri_dist(0,9);
	char str;
	for(int count=1;count<=20;count++){
		str = 'a'+prior_dist(mt);
		f<<str;
	}
	for(int count=1;count<=5;count++){
		str = '0'+pri_dist(mt);
		f<<str;
	}

	cout<<endl;
	//cin.sync(); cin.get();
}

int main(){
	password();
	return 0;
}