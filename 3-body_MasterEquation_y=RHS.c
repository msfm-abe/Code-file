#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double factorial(int k);
double RightHandSide(double prob,double c);

#define TEMPERATURE 0.5

int main(){

	double prob;
	double c;
	char filename[100];

	FILE *fp;
	
	for(c=3.00;c<=5.01;c+=0.01){
		sprintf(filename,"3-body_MasterEquation_y=RHS_T=0.5_c=%d_div100.dat",(int)(c*100));

		fp = fopen(filename,"wb");
		for(prob=0.0;prob<=1.01;prob+=0.01){
			fprintf(fp,"%f %f\n",prob,RightHandSide(prob,c));
		}

		fclose(fp);
	}

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

	return -prob+summation_l;
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