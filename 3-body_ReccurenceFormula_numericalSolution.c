#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define EPS pow(10.0,-8.0)
#define NMAX 100
#define TEMPERATURE 0.5

double factorial(int k);
double RightHandSide(double prob,double c);
double dif_RightHandSide(double prob,double c);

int main(){

	int n;
	double prob;
	double c;
	double initial_prob;
	double d;

	/*FILE *fp;
	fp = fopen("test8_2_sono2.dat","wb");*/

	//c = 3.76;

	for(initial_prob=1.0;initial_prob>=0.1;initial_prob-=0.1){
		prob = initial_prob;
		n = 0;

		for(c=3.0;c<=5.0;c+=0.1){

			do{
				d = -RightHandSide(prob,c)/dif_RightHandSide(prob,c);
				prob = prob+d;
				n++;
			}while(fabs(d)>EPS && n<NMAX);

			if(n==NMAX){
				printf("fault\n");
			}else{
				printf("initialProb %f c=%f %f\n",initial_prob,c,prob);
			}

		}
		printf("\n");

	}

	//fclose(fp);

	getchar();

	return 0;
}

//double RightHandSide(double m,double T){
//	return m-tanh(m/T);
//}
//
//double dif_RightHandSide(double m,double T){
//	return 1.0-(1.0/pow(cosh(m/T),2.0))*(1/T);
//}

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

double dif_RightHandSide(double prob,double c){

	int k,l;
	double summation_k;
	double summation_l;

	summation_l = 0.0;

	for(l=0;l<=20;l++){

		summation_k = 0.0;

		for(k=0;k<=l;k++){
			summation_k += ( factorial(l)/( factorial(k)*factorial(l-k) ) )
				*( -2.0*(2.0*prob-1.0)*pow(2.0*pow(prob,2)-2.0*prob+1.0,k-1)*pow(2.0*prob-2.0*pow(prob,2),l-k-1)*(2.0*prob*(prob-1.0)*l+l-k)
				/( 1.0+exp(-2.0*(2.0*k-l)/TEMPERATURE) ) );
		}

		summation_l += ( exp(-c)*pow(c,l)/factorial(l) )*summation_k;
	}

	return -1.0+summation_l;

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