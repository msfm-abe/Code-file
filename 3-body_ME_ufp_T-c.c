#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define mag 1.0
#define EPS pow(10.0,-8.0)
#define NMAX 10

double RightHandSide(double temperature,double c);
double dRightHandSide(double temperature,double c);
double factorial(int k);

int main(){

	int n;
	double temperature;
	double d;
	double c;
	double temperaturemin,temperaturemax;
	double y1,y2;
	double h;
	int bunkatsusu;

	FILE *fp;
	fp = fopen("unstable.dat","wb");

	for(c=0.0;c<=10.0;c+=0.01){

		temperaturemin = 0.01;
		temperaturemax = 1.01;
		bunkatsusu = 100;

		h = (temperaturemax-temperaturemin)/bunkatsusu; // h=0.01

		y1 = RightHandSide(temperaturemin,c);

		for(temperature=temperaturemin+h;temperature<=temperaturemax;temperature+=h){

			if(mag>=1.0){	
				fprintf(fp,"%f 0.0\n",c);
				break;
			}

			y2 = RightHandSide(mag,c);

			if(y1*y2<0){
				n = 0;
				do{
					d = -RightHandSide(temperature,c)/dRightHandSide(temperature,c);
					temperature = temperature+d;
					n++;
				}while(fabs(d)>EPS && n<NMAX);

				if(n==NMAX){
					fprintf(fp,"%f fault\n",c);
				}else{
					fprintf(fp,"%f %f\n",c,fabs(temperature));
				}
				break;
			}else{
				y1 = y2;
			}

		}//mag

	}//c

	fclose(fp);
	//getchar();
	return 0;

}

double RightHandSide(double temperature,double c){
 
    int k,l;
    double summation_k;
    double summation_l;
 
    summation_l = 0.0;
 
    for(l=0;l<=30;l++){
 
        summation_k = 0.0;
 
        for(k=0;k<=l;k++){
            summation_k += ( factorial(l)/( factorial(k)*factorial(l-k) ) )
                *( pow(pow((1+mag)/2.0,2)+pow((1-mag)/2.0,2),k)*pow(0.5*(1+mag)*(1-mag),l-k)/( 1.0+exp(-2.0*(2.0*k-l)/temperature) ) );

		}
 
        summation_l += ( exp(-c)*pow(c,l)/factorial(l) )*summation_k;
    }
 
    return -(mag+1.0)+2.0*summation_l;
 
}

double dRightHandSide(double temperature,double c){
	
	int k,l;
	double summation_k,summation_l;

	summation_l = 0.0;
	
	for(l=0;l<=30;l++){
		summation_k = 0.0;
		for(k=0;k<=l;k++){
			summation_k += ( factorial(l)/(factorial(k)*factorial(l-k)) )*
				( pow( pow((1+mag)/2.0,2)+pow((1-mag)/2.0,2),k-1 ) )*( pow((1-pow(mag,2))/2.0,l-k-1) )*( k-l*(1+pow(mag,2))/2.0 )
				/( 1+exp(-2*(2*k-l)/TEMPERATURE) );
		}
		summation_l += ( exp(-c)*pow(c,l)/factorial(l) )*summation_k;
	}

	return -1.0+2.0*mag*summation_l;
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