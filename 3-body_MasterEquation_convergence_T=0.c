#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define c 5.0

double RightHandSide(double probability);
double factorial(int k);

int main(){

	int i,j,k;
	double initialCondition;
	double probability;
	double time;
	double h,a,b,n;
	double k1,k2,k3,k4;

	FILE *fp;
	fp = fopen("3-body_InitialCondition_Convergence_T=0_c=5.0_Correction.dat","wb");

	//for(i=5;i>=0;i--)
	//{
	//	factorial(i);
	//}

	b = 300.0;
	a = 0.0;
	n = 3000.0;
	
	h = (b-a)/n; //h=0.1‚Éfix
	//printf("%f\n",h);
	
	for(initialCondition=0.0;initialCondition<=1.0;initialCondition+=0.1)
	{
		probability = initialCondition;
		for(time=0;time<=b;time+=h)
		{
			//printf("%f %f\n",time,probability);
			fprintf(fp,"%f %f\n",time,probability);

			k1 = RightHandSide(probability);
			k2 = RightHandSide(probability+h*k1/2);
			k3 = RightHandSide(probability+h*k2/2);
			k4 = RightHandSide(probability+h*k3);

			probability = probability+h*(k1+k2+k3+k4)/6.0;

			//time += h;
		}
		//printf("\n");
		fprintf(fp,"\n");
	}

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

double RightHandSide(double probability){

	double summation_l=0,summation_k=0;
	int k,l;

	summation_l = 0.0;
	for(l=0;l<=20;l++)
	{
		summation_k = 0.0;
		for(k=0;k<=l;k++)
		{
			if((2*k-l)>0){
				summation_k += ( factorial(l)/( factorial(k)*factorial(l-k) ) )*( pow(2*probability*probability-2*probability+1,k)*pow(2*probability-2*probability*probability,l-k) );
			}else if((2*k-l)==0){
				summation_k += 0.5*( factorial(l)/( factorial(k)*factorial(l-k) ) )*( pow(2*probability*probability-2*probability+1,k)*pow(2*probability-2*probability*probability,l-k) );
			}else if((2*k-l)<0){
				summation_k += 0.0;
			}
  			//summation_k += ( factorial(l)/( factorial(k)*factorial(l-k) ) )*( pow(2*probability*probability-2*probability+1,k)*pow(2*probability-2*probability*probability,l-k) )/( 1.0+exp(-4*(k-l/2.0)/temperature) );
			//summation_k += ( factorial(l)/( factorial(k)*factorial(l-k) ) )*( pow(probability,k)*pow(1-probability,l-k) )/( 1.0+exp(-4*(k-l/2.0)/temperature) );
		}
		summation_l += ( exp(-c)*pow(c,l)/factorial(l) )*summation_k;
	}

	return (-probability+summation_l);
}