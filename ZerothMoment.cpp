#include<iostream>
#include<math.h>
#include<stdio.h>
#include<cstdlib>

using namespace std;

static const int Row = 999;
static const double delta_x = 0.003;

double CalcMoment0(FILE *fp){

	double moment0 = 0.0;
	double dataPoint[Row]={0.0};
	double time;
	int count=0;

	for(int i=0;i<2*Row;i++){
		if(i%2==0){
			fscanf(fp,"%lf\t",&time);
		}else{
			fscanf(fp,"%lf",&dataPoint[count]);
			//cout<<"dataPoint["<<count<<"]="<<dataPoint[count]<<endl;
			count++;
		}
	}

	for(int i=0;i<Row-1;i++)
		moment0 += delta_x*((dataPoint[i+1]+dataPoint[i])/2.0);

	cout<<"M0="<<moment0<<endl;

	return moment0;

}

int main(){

	FILE *fp1,*fp2,*fp3,*fp4;
	fp1 = fopen("Eq21EF_alpha=1.0e-03_N=100_T=3.0_width=11.dat","rb"); if(fp1==NULL) exit(1);
	fp2 = fopen("Eq21EF_alpha=1.0e-03_N=100_T=3.0_width=21.dat","rb"); if(fp2==NULL) exit(1);
	fp3 = fopen("Eq21EF_alpha=1.0e-03_N=100_T=3.0_width=31.dat","rb"); if(fp3==NULL) exit(1);
	fp4 = fopen("StandardDiffusionCurve.dat","rb"); if(fp4==NULL) exit(1);

	cout<<"sdc ";CalcMoment0(fp4);
	cout<<"width=10 ";CalcMoment0(fp1);
	cout<<"width=20 ";CalcMoment0(fp2);
	cout<<"width=30 ";CalcMoment0(fp3);

	fclose(fp1);fclose(fp2);fclose(fp3);
	cin.sync(); cin.get();
	return 0;
}