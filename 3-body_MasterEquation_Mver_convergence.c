#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define C 3.9
#define TEMPERATURE 0.5

double RightHandSide(double probability,double temperature,double c);
double factorial(int k);

int main(){

	int i,j,k;
	double initialCondition;
	double probability;
	double time;
	double h,a,b,n;
	double k1,k2,k3,k4;
	double temperature = TEMPERATURE;
	double c = C;
	char filename[100];

	FILE *fp;

	fp = fopen("3-body_MasterEq_M.ver_TimeEvolution_T=0.5_c=3.9.dat","wb");


	b = 1000.0;
	a = 0.0;
	n = 10000.0;
	
	h = (b-a)/n;//h=0.1‚Éfix
	
	for(initialCondition=0.0;initialCondition<=1.0;initialCondition+=0.1){

		probability = initialCondition;

		for(time=0;time<=b;time+=h)
		{
			fprintf(fp,"%f %f\n",time,2*probability-1);

			k1 = RightHandSide(probability,temperature,c);
			k2 = RightHandSide(probability+h*k1/2,temperature,c);
			k3 = RightHandSide(probability+h*k2/2,temperature,c);
			k4 = RightHandSide(probability+h*k3,temperature,c);

			probability = probability+h*(k1+2*k2+2*k3+k4)/6.0;

		}

		fprintf(fp,"\n");

	}//initial condition

	fclose(fp);

	//getchar();
	return 0;
}

double RightHandSide(double probability,double temperature,double c){

	double summation_l=0,summation_k=0;
	int k,l;

	summation_l = 0.0;
	for(l=0;l<=20;l++)
	{
		summation_k = 0.0;
		for(k=0;k<=l;k++)
		{

			summation_k += ( factorial(l)/(factorial(k)*factorial(l-k)) )*( pow(pow(probability,2)+pow(1-probability,2),k)
				           *pow(2*probability*(1-probability),l-k) )/( 1+exp(-2*(2*k-l)/temperature) );

		}
		summation_l += ( exp(-c)*pow(c,l)/factorial(l) )*summation_k;
	}

	return (-probability+summation_l);
}

double factorial(int k){

	double Factorial;
	int i;

	Factorial = 1;

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