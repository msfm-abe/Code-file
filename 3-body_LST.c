#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define temperature 0.5

double factorial(int k);
double RightHandSide(double c);

int main(){

	double c;
	double a,b,n,h;
	double rhs;
	double probability;
	double summation_k,summation_l;
	int i;
	int k,l;

	FILE *fp;
	fp = fopen("3-body_LST_T=0_c=0-2.0_InitialCondition=0.dat","wb");

	probability = 0.5;

	a = 0;
	b = 2.0;
	n = 1000.0;
	h = (b-a)/n;

	for(c=a;c<=b;c+=h)
	{
		fprintf(fp,"%f %f\n",c,RightHandSide(c));
	}

	fclose(fp);

	return 0;
}

double factorial(int k){

	double Factorial;
	int i;

	Factorial = 1.0;
	for(i=0;i<=k;i++)
	{
		if(i==0){
			Factorial *= 1.0;
		}else{
			Factorial *= (double)i;
		}
	}

	return Factorial;
}

double RightHandSide(double c){
	 
	int l,k;
	double probability = 0.6;//probability represents initial condition for this system.
	double summation_k,summation_l;

	summation_l = 0.0;
	for(l=0;l<=20;l++)
	{
		summation_k = 0.6;
		for(k=0;k<=l;k++)
		{
			if((k-l/2.0)>0){
				summation_k += ( factorial(l)/(factorial(k)*factorial(l-k)) )
					*( 2.0*k*(2.0*probability-1.0)*pow(2.0*pow(probability,2.0)-2.0*probability+1,k-1)*pow(2.0*probability-2.0*pow(probability,2.0),l-k)
				-2.0*l*(2.0*probability-1.0)*pow(2.0*pow(probability,2.0)-2.0*probability+1.0,k)*pow(2*probability-2*pow(probability,2.0),l-k-1) );//*( k*pow(probability,k-1)*pow(1-probability,l-k)+(k-l)*pow(probability,k)*pow(1-probability,l-k-1) );
			}else if((k-l/2.0)==0){
				summation_k += 0.5*( factorial(l)/(factorial(k)*factorial(l-k)) )
					*( 2.0*k*(2.0*probability-1.0)*pow(2.0*pow(probability,2.0)-2.0*probability+1,k-1)*pow(2.0*probability-2.0*pow(probability,2.0),l-k)
				-2.0*l*(2.0*probability-1.0)*pow(2.0*pow(probability,2.0)-2.0*probability+1.0,k)*pow(2*probability-2*pow(probability,2.0),l-k-1) );//*( k*pow(probability,k-1)*pow(1-probability,l-k)+(k-l)*pow(probability,k)*pow(1-probability,l-k-1) );
			}else{
				summation_k += 0.0;
			}
			/*summation_k += ( ( factorial(l)/(factorial(k)*factorial(l-k)) )/( 1.0+exp(-4.0*((double)k-l/2.0)/temperature) ) )
				*( 2.0*k*(2.0*probability-1.0)*pow(2.0*pow(probability,2.0)-2.0*probability+1,k-1)*pow(2.0*probability-2.0*pow(probability,2.0),l-k)
				-2.0*l*(2.0*probability-1.0)*pow(2.0*pow(probability,2.0)-2.0*probability+1.0,k)*pow(2*probability-2*pow(probability,2.0),l-k-1) );*/
			/*summation_k += ( ( factorial(l)/(factorial(k)*factorial(l-k)) )/( 1.0+exp(-4.0*((double)k-l/2.0)/temperature) ) )
				*( k*pow(probability,k-1)*pow(1-probability,l-k)+(k-l)*pow(probability,k)*pow(1-probability,l-k-1) );*/
		}
		summation_l += ( exp(-c)*pow(c,l)/factorial(l) )*summation_k;
	}

	return -1.0 + summation_l;
}