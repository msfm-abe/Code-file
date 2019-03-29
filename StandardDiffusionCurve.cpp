#include<iostream>
#define _USE_MATH_DEFINES
#include<math.h>
#include<stdio.h>
#include<cstdlib>
#include<fstream>
using namespace std;

void SDC(){
	double flow;
	ofstream f("StandardDiffusionCurve.dat");

	for(double time=0.003;time<3;time+=0.003){
		
		flow = 0.0;
		for(int n=0;n<=1000;n++)
			flow += M_PI*pow(-1.0,n)*(2.0*n+1)*exp(-pow(n+0.5,2)*pow(M_PI,2)*time);

		f<<time<<" "<<flow<<endl;
	}
}

int main(){
	SDC();
	return 0;
}