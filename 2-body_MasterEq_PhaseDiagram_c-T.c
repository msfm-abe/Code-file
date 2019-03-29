#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double RightHandSide(double mag,double c,double temperature);
double factorial(int k);

int main(){
	
	int i,j,k;
	double InitialMagnetization;
	double mag;
	double temperature;
	double c;
	double time;
	double h,a,b;
	int n;
	double k1,k2,k3,k4;

	FILE *fp;

	fp = fopen("2-body_c-T_PhaseDiagram.dat","wb");

	b = 1000.0;
	a = 0.0;
	n = 10000;
	
	h = (b-a)/n; //h=0.1

	InitialMagnetization = 1.0;

	for(temperature=0.0;temperature<=2.01;temperature+=0.01){

		for(c=1.500;c<=3.000;c+=0.001){

			mag = InitialMagnetization;

			for(time=a;time<b;time+=h){

				k1 = RightHandSide(mag,c,temperature);
				k2 = RightHandSide(mag+k1*h/2,c,temperature);
				k3 = RightHandSide(mag+k2*h/2,c,temperature);
				k4 = RightHandSide(mag+k3*h,c,temperature);

				mag = mag+h*(k1+2*k2+2*k3+k4)/6;

			}

			if(mag>0.0){
				fprintf(fp,"%f %f\n",c,temperature);
				break;
			}
	
		}//c
	
	}//temperature

	fclose(fp);
			
	return 0;
}

double RightHandSide(double mag,double c,double temperature){
	
	double summation_k;
	double summation_l;
	int k,l;

	summation_l = 0.0;

	for(l=0;l<=20;l++){

		summation_k = 0.0;

		for(k=0;k<=l;k++){

			if(temperature>0.0){
				summation_k += ( factorial(l)/(factorial(k)*factorial(l-k)) )*
					( pow((1+mag)/2,k)*pow((1-mag)/2,l-k) )/( 1+exp(-2*(2*k-l)/temperature) );
			}else{
				if(2*k-l>0){
					summation_k += ( factorial(l)/(factorial(k)*factorial(l-k)) )*
					( pow((1+mag)/2,k)*pow((1-mag)/2,l-k) );
				}else if(2*k-l==0){
					summation_k += 0.5*( factorial(l)/(factorial(k)*factorial(l-k)) )*
					( pow((1+mag)/2,k)*pow((1-mag)/2,l-k) );
				}else if(2*k-l<0){
					summation_k += 0.0;
				}
			}

		}//k

		summation_l += ( exp(-c)*pow(c,l)/factorial(l) )*summation_k;

	}//l

	return -(1+mag)+2*summation_l;
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