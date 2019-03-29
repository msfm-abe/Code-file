#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define TEMPERATURE 0.5

double factorial(int k);
double RHS(double mag,double c);

int main(){

	double mag;
	double a,b,n,h;
	double c;
	char filename[100];

	FILE *fp;

	a = 0.0;
	b = 1.01;
	n = 101.0;
	h = (b-a)/n;


	for(c=3.00;c<4.51;c+=0.01){

		sprintf(filename,"3-body_MasterEquation_M.ver_y=m_y=f_T-0.5_c=%d_div100.dat",(int)(c*100));
		fp = fopen(filename,"wb");

		for(mag=a;mag<=b;mag+=h){
			fprintf(fp,"%f %f\n",mag,RHS(mag,c));
		}

		fclose(fp);
	}

	return 0;
}

double factorial(int k){

	int i;
	double product;

	product = 1.0;

	for(i=k;i>=0;i--){
		if(i!=0){
			product *= (double)i;
		}else if(i==0){
			product *= 1.0;
		}
	}

	return product;
}

double RHS(double mag,double c){

	int k,l;
	double summation_k,summation_l;

	summation_l = 0.0;
	for(l=0;l<=20;l++){
		summation_k = 0.0;
		for(k=0;k<=l;k++){
			summation_k += ( factorial(l)/(factorial(k)*factorial(l-k)) )
				*( pow(pow((1+mag)/2.0,2)+pow((1-mag)/2.0,2),k)*pow(0.5*(1+mag)*(1-mag),l-k)/( 1.0+exp(-2.0*(2.0*k-l)/TEMPERATURE) ) );
		}
		summation_l += ( pow(c,l)*exp(-c)/factorial(l) )*summation_k;
	}

	return -1.0+2.0*summation_l;
}