#include<iostream>
#define _USE_MATH_DEFINES
#include<stdio.h>
#include<math.h>
#include<cstdlib>
#include<fstream>
using namespace std;

static const double L = 1.0; //定義域
static const int N = 100; //空間方向の分割数
static const double Time = 4.0; //計算時間
static const double Diff = 1.0; //拡散定数

double* InitialCondition(double *u,double delta_x,int width);

int main(){

	double alpha = 0.001; //(delta_t*Diff)/(delta_x)^2
	double delta_x = L/N;
	double delta_t = 0;
	double u[N+3]={0.0};
	char name[100];
	int count = 0;
	int width = 31; //初期条件に入れるデルタ関数の幅

	delta_t = alpha*pow(delta_x,2)/Diff;

	int M = Time/delta_t; cout<<M<<endl; //時間方向の分割数を計算
//	cout<<"M*delta_t="<<M*delta_t<<endl;

	InitialCondition(u,delta_x,width);

	cout<<"delta_t = "<<delta_t<<endl;
	cout<<"delta_x = "<<delta_x<<" (D*delta_t)/pow(delta_x,2)="<<(Diff*delta_t)/pow(delta_x,2)<<endl;

	double t = 0.0;

	sprintf(name,"eq21EF_alpht=%2.1e_Diff=%3.1f_delta_x=%3.1e_delta_t=%3.1e_width=%d.dat",alpha,Diff,delta_x,delta_t,width-1);
	ofstream f(name);

	while(t<=Time){
		//if(count==(int)(M*0.1) || count==M/2 || count==M){
		//	sprintf(name,"eq14_alpht=%3.1e_Diff=%3.1f_delta_x=%3.1e_delta_t=%3.1e_time=%f_width=%d.dat",alpha,Diff,delta_x,delta_t,t,width-1);
		//	ofstream f(name);
		//	for(int j=1;(j-0.5)*delta_x<=L;j++) f<<(j-0.5)*delta_x<<" "<<t<<" "<<u[j]<<endl;
		//		cout<<"yeah"<<endl;
		//}

		for(double k=0.01;k<=1;k+=0.01){
			if(count==(M*k)){
				f<<t<<" "<<-(u[N+1]-u[N-1])/delta_x<<endl;
			}
		}

		for(int i=1;i<=N+1;i++)
			u[i] = u[i] + ((Diff*delta_t)/pow(delta_x,2))*(u[i-1]-2.0*u[i]+u[i+1]);
		u[0] = u[1]; u[N+2] = -u[N+1];

		//u[0] = -u[1]; u[N+2] = -u[N+1];

		t += delta_t;
		count++;
	}

	printf("finish");

	//cin.sync(); cin.get();
	return 0;
}

double* InitialCondition(double *u,double delta_x,int width){

	for(int i=1;i<=N+1;i++){
		if(i<=width)
			u[i] = 1.0/(delta_x*(width-1));
		else u[i] = 0.0;
	}
	u[0] = u[1]; u[N+2] = -u[N+1];

	//for(int i=1;i<=N+1;i++) 
	//	*(u+i) = 2.0*sin((i-0.5)*delta_x)+sin(2.0*(i-0.5)*delta_x);

	//u[0] = -u[1]; u[N+2] = -u[N+1];
	
	return u;
}