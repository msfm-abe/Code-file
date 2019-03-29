#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define temperature 0.5

double RightHandSide(double probability,double c);
double factorial(int k);

int main(void){

	int i,j,k;
	double probability;
	double k1,k2,k3,k4;
	double h,a,b,n;
	double time;
	double c;

	FILE *fp;
	fp = fopen("3-body_MasterEquation_c_izon_T=0.5_InitialProbability=0.95.dat","wb");


	a = 0.0;
	b = 300.0;
	n = 3000.0;
	h = (b-a)/n;//h is fixed to 0.1;


	for(c=3.0;c<=5.0;c+=0.01){

		probability = 0.95;

		for(time=a;time<=b;time+=h){
			k1 = RightHandSide(probability,c);
			k2 = RightHandSide(probability+h*k1/2.0,c);
			k3 = RightHandSide(probability+h*k2/2.0,c);
			k4 = RightHandSide(probability+h*k3,c);

			probability = probability + h*(k1+2.0*k2+2.0*k3+k4)/6.0;
		}

		fprintf(fp,"%f %f\n",c,2.0*probability-1.0);
	}

	fclose(fp);

	return 0;
}

double RightHandSide(double probability,double c){

	int k,l;
	double summation_k;
	double summation_l;

	summation_l = 0.0;

	for(l=0;l<=20;l++){

		summation_k = 0.0;

		for(k=0;k<=l;k++){
			summation_k += ( factorial(l)/( factorial(k)*factorial(l-k) ) )
				*( pow(2.0*pow(probability,2)-2.0*probability+1.0,k)*pow(2.0*probability-2.0*pow(probability,2),l-k)/( 1.0+exp(-2.0*(2.0*k-l)/temperature) ) );
		}

		summation_l += ( exp(-c)*pow(c,l)/factorial(l) )*summation_k;
	}

	return -probability+summation_l;

}

double factorial(int k){

	int i;
	double product=1.0;

	for(i=k;i>=0;i--){
		if(i!=0){
			product *= (double)i;
		}else if(i==0){
			product *= 1.0;
		}
	}

	return product;
}