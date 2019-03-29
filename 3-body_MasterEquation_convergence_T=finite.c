#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//#define c 3.75
#define temperature 0.5

double RightHandSide(double probability,double c);
double factorial(int k);

int main(){

	int i,j,k;
	double initialCondition;
	double probability;
	double time;
	double h,a,b,n;
	double k1,k2,k3,k4;
	double c;
	char filename[100];

	FILE *fp;

	for(c=3.7;c<=3.9;c+=0.1){
		sprintf(filename,"3-body_InitialCondition=0.95_Convergence_T=0.5_c=%f.dat",c);

	fp = fopen(filename,"wb");

	//for(i=5;i>=0;i--)
	//{
	//	factorial(i);
	//}

	b = 300.0;
	a = 0.0;
	n = 3000.0;
	
	h = (b-a)/n;//h=0.1‚Éfix
	//printf("%f\n",h);
	
	//for(initialCondition=0.0;initialCondition<=1.0;initialCondition+=0.1)
	//{
		//probability = initialCondition;
	probability = 0.95;
		for(time=0;time<=b;time+=h)
		{
			//printf("%f %f\n",time,probability);
			fprintf(fp,"%f %f\n",time,probability);

			k1 = RightHandSide(probability,c);
			k2 = RightHandSide(probability+h*k1/2,c);
			k3 = RightHandSide(probability+h*k2/2,c);
			k4 = RightHandSide(probability+h*k3,c);

			probability = probability+h*(k1+2*k2+2*k3+k4)/6.0;

			//time += h;
		}
	}
		//printf("\n");
		fprintf(fp,"\n");
	//}

	fclose(fp);

	//getchar();
	return 0;
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

	//printf("%f\n",Factorial);

	return Factorial;
}

double RightHandSide(double probability,double c){

	double summation_l=0,summation_k=0;
	int k,l;

	summation_l = 0.0;
	for(l=0;l<=20;l++)
	{
		summation_k = 0.0;
		for(k=0;k<=l;k++)
		{
			//summation_k += ( factorial(l)/( factorial(k)*factorial(l-k) ) )*( pow(pow(probability,3)+3*probability*pow(1-probability,2),k)
				//*pow(pow(1-probability,3)+3*pow(probability,2)*(1-probability),l-k) )/( 1.0+exp(-4*(k-l/2.0)/temperature) );
  			summation_k += ( factorial(l)/( factorial(k)*factorial(l-k) ) )*( pow(2*probability*probability-2*probability+1,k)*pow(2*probability-2*probability*probability,l-k) )/( 1.0+exp(-4*(k-l/2.0)/temperature) );
			//summation_k += ( factorial(l)/( factorial(k)*factorial(l-k) ) )*( pow(probability,k)*pow(1-probability,l-k) )/( 1.0+exp(-4*(k-l/2.0)/temperature) );
		}
		summation_l += ( exp(-c)*pow(c,l)/factorial(l) )*summation_k;
	}

	return (-probability+summation_l);
}