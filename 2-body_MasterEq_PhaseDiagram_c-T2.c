#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double rhs(double probability,double c,double temperature);
double drhs(double probability,double c,double temperature);
double factorial(int k);

int main(void){

	double probability;
	double temperature;
	double c;
	double epsilon=0.0;
	FILE *fp;
	fp = fopen("testsouzu.dat","wb");

	for(temperature=0.0;temperature<=1.6;temperature+=0.001){

		for(c=1.7;c<=2.6;c+=0.001){

			if(drhs(0.5,c,temperature)>=0.0){
				fprintf(fp,"%f %f\n",c,temperature);
				break;
			}
		}
	}

	fclose(fp);
	//getchar();
	return 0;
}

double rhs(double probability,double c,double temperature){

	int k,l;
	double summation_k,summation_l;

	summation_l = 0.0;

	for(l=0;l<=20;l++){

		summation_k = 0.0;

		for(k=0;k<=l;k++){

			if(temperature!=0.0){
				summation_k += ( factorial(l)/(factorial(k)*factorial(l-k)) )
					*( pow(probability,k)*pow(1.0-probability,l-k) )/( 1.0+exp(-2*(2.0*k-l)/temperature) );
			}else if(temperature==0.0){
				if((2*k-l)>0){
					summation_k += ( factorial(l)/(factorial(k)*factorial(l-k)) )
						*( pow(probability,k)*pow(1.0-probability,l-k) );
				}else if((2*k-l)==0){
					summation_k += 0.5*( factorial(l)/(factorial(k)*factorial(l-k)) )
						*( pow(probability,k)*pow(1.0-probability,l-k) );
				}else if((2*k-l)<0){
					summation_k += 0.0;
				}
			}

		}//k

		summation_l += ( exp(-c)*pow(c,l)/factorial(l) )*summation_k;

	}//l

	return -probability+summation_l;
}

double drhs(double probability,double c,double temperature){

	int k,l;
	double summation_k,summation_l;

	summation_l = 0.0;

	for(l=0;l<=20;l++){

		summation_k = 0.0;

		for(k=0;k<=l;k++){

			if(temperature!=0.0){
				summation_k += ( factorial(l)/(factorial(k)*factorial(l-k)) )
					*( (k-l*probability)*pow(probability,k-1)*pow(1.0-probability,l-k-1) )/( 1.0+exp(-2*(2.0*k-l)/temperature) );
			}else if(temperature==0.0){
				if((2*k-l)>0){
					summation_k += ( factorial(l)/(factorial(k)*factorial(l-k)) )
						*( (k-l*probability)*pow(probability,k)*pow(1.0-probability,l-k) );
				}else if((2*k-l)==0){
					summation_k += 0.5*( factorial(l)/(factorial(k)*factorial(l-k)) )
						*( (k-l*probability)*pow(probability,k)*pow(1.0-probability,l-k) );
				}else if((2*k-l)<0){
					summation_k += 0.0;
				}
			}

		}//k

		summation_l += ( exp(-c)*pow(c,l)/factorial(l) )*summation_k;

	}//l

	return -1.0+summation_l;
}

double factorial(int k){

	int i,j;
	double product;

	product = 1.0;

	for(i=0;i<=k;i++){

		product *= (i!=0)? (double)i:1.0;

	}

	return product;
}