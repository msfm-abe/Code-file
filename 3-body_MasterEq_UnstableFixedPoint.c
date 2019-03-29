#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define TEMPERATURE 0.5
#define EPS pow(10.0,-8.0)
#define NMAX 10

double RightHandSide(double mag,double c);
double dRightHandSide(double mag,double c);
double factorial(int k);

int main(){

	int n;
	double mag;
	double d;
	double c;
	double magmin,magmax;
	double y1,y2;
	double h;
	int bunkatsusu;

	FILE *fp;
	fp = fopen("unstable.dat","wb");

	for(c=0.0;c<=10.0;c+=0.001){

		magmin = 0.01;
		magmax = 1.01;
		bunkatsusu = 1000;

		h = (magmax-magmin)/bunkatsusu; // h=0.01

		y1 = RightHandSide(magmin,c);

		for(mag=magmin+h;mag<=magmax;mag+=h){

			if(mag>=1.0){	
				fprintf(fp,"%f 0.0\n",c);
				break;
			}

			y2 = RightHandSide(mag,c);

			if(y1*y2<0){
				n = 0;
				do{
					d = -RightHandSide(mag,c)/dRightHandSide(mag,c);
					mag = mag+d;
					n++;
				}while(fabs(d)>EPS && n<NMAX);

				if(n==NMAX){
					fprintf(fp,"%f fault\n",c);
				}else{
					fprintf(fp,"%f %f\n",c,fabs(mag));
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

double RightHandSide(double mag,double c){
 
    int k,l;
    double summation_k;
    double summation_l;
 
    summation_l = 0.0;
 
    for(l=0;l<=30;l++){
 
        summation_k = 0.0;
 
        for(k=0;k<=l;k++){
            summation_k += ( factorial(l)/( factorial(k)*factorial(l-k) ) )
                *( pow(pow((1+mag)/2.0,2)+pow((1-mag)/2.0,2),k)*pow(0.5*(1+mag)*(1-mag),l-k)/( 1.0+exp(-2.0*(2.0*k-l)/TEMPERATURE) ) );

		}
 
        summation_l += ( exp(-c)*pow(c,l)/factorial(l) )*summation_k;
    }
 
    return -(mag+1.0)+2.0*summation_l;
 
}

double dRightHandSide(double mag,double c){
	
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