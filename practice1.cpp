#include<iostream>
#define _USE_MATH_DEFINES
#include<stdio.h>
#include<math.h>
#include<cstdlib>
#include<fstream>
using namespace std;

static const int M = 10000000;
static const double T = 1.0;
static const double L = M_PI;
static const double N = M_PI*100;
static const int intN = (int)M_PI*100;
static const double D = 1.0;

double* InitialCondition(double *u,double delta_x);

int main(){

	double delta_x = L/intN;
	double delta_t = T/M;
	double u[intN+2]={0.0};
	char name[30];
	int count = 0;

	InitialCondition(u,delta_x);

	cout<<"delta_t = "<<delta_t<<endl;
	cout<<"delta_x = "<<delta_x<<" (D*delta_t)/pow(delta_x,2)="<<(D*delta_t)/pow(delta_x,2)<<endl;

	double t = 0.0;

	if( (D*delta_t)/pow(delta_x,2)<0.5 ){
		printf("yes\n");
		for(count=0;count<=M;count++){
			if(count==0 || count==1*(M/10) || count==2*(M/10) || count==3*(M/10) || count==4*(M/10) || count==5*(M/10) || count==6*(M/10) ||
				count==7*(M/10) || count==8*(M/10) || count==9*(M/10) ||count==10*(M/10)){
					sprintf(name,"practice1_time=%f_.dat",t);
					ofstream f(name);
					for(int j=1;j<=intN;j++) f<<(j-0.5)*delta_x<<" "<<t<<" "<<u[j]<<endl;
			}

			for(int i=1;i<=intN;i++)
				u[i] = u[i] + ((D*delta_t)/pow(delta_x,2))*(u[i-1]-2.0*u[i]+u[i+1]);
			u[0] = -u[1]; u[intN+1] = -u[intN];

			t += delta_t;

		}
		printf("finish");
	}else{
		cout<<"No"<<endl;
	}

	cin.sync(); cin.get();
	return 0;
}

double* InitialCondition(double *u,double delta_x){
	for(int i=1;i<=intN;i++) *(u+i) = 2.0*sin((i-0.5)*delta_x)+sin(2.0*(i-0.5)*delta_x);
	u[0] = -u[1]; u[intN+1] = -u[intN];
	return u;
}