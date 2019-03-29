#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define TEMPERATURE 0.5

double factorial(int k);
double RightHandSide(double prob,double c);

int main(){

	int n;
	double prob;
	double c;
	double h,a,b;
	char filename[100];

	FILE *fp;

	a = 0.0;
	b = 1.01;
	n = 101;
	h = (b-a)/n; //h is fixed to 0.01

	for(c=3.00;c<=4.0;c+=0.1){
		sprintf(filename,"3-body_MasterEquation_NumericalSolution_c=%d_div100.dat",(int)(c*100));
		fp = fopen(filename,"wb");

		for(prob=a;prob<=b;prob+=h){
			fprintf(fp,"%f %f\n",prob,RightHandSide(prob,c));
		}

		fclose(fp);
	}

	//getchar();

	return 0;
}

double RightHandSide(double prob,double c){

	int k,l;
	double summation_k;
	double summation_l;

	summation_l = 0.0;

	for(l=0;l<=20;l++){

		summation_k = 0.0;

		for(k=0;k<=l;k++){
			summation_k += ( factorial(l)/( factorial(k)*factorial(l-k) ) )
				*( pow(2.0*pow(prob,2)-2.0*prob+1.0,k)*pow(2.0*prob-2.0*pow(prob,2),l-k)/( 1.0+exp(-2.0*(2.0*k-l)/TEMPERATURE) ) );
		}

		summation_l += ( exp(-c)*pow(c,l)/factorial(l) )*summation_k;
	}

	return summation_l;

}

double factorial(int k){
	double Factorial;
	int i;

	Factorial = 1.0;

	for(i=k;i>=0;i--)
	{
		if(i!=0){
			Factorial *= (double)i;
		}else if(i==0){
			Factorial *= (double)1;
		}
	}

	return Factorial;
}