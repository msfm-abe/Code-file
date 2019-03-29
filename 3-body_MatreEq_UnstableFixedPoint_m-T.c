#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define C 4.0
#define EPS pow(10.0,-8.0)
#define NMAX 1000

double RightHandSide(double mag,double temperature,double c);
double dRightHandSide(double mag,double temperature,double c);
double factorial(int k);

int main(){

	int n;
	double mag;
	double d;
	double temperature;
	double c=C;
	double magmin,magmax;
	double y1,y2;
	double h;
	int bunkatsusu;
	int hantei;
	char filename[100];

	FILE *fp;
	fp = fopen("3-body_MasterEq_UnstableFixedPoint_m-T_c=10.0.dat","wb");

	for(temperature=0.0;temperature<=2.01;temperature+=0.01){

		magmin = 0.01;
		magmax = 1.01;
		bunkatsusu = 1000; //h=0.001

		h = (magmax-magmin)/bunkatsusu; // h=0.01

		y1 = RightHandSide(magmin,temperature,c);

		for(mag=magmin+h;mag<=1.01;mag+=h){

			if(mag>=1.0){	
				//fprintf(fp,"%f 0.0\n",temperature);
				break;
			}

			y2 = RightHandSide(mag,temperature,c);

			if(y1*y2<0){
				n = 0;
				do{
					d = -RightHandSide(mag,temperature,c)/dRightHandSide(mag,temperature,c);
					mag = mag+d;
					n++;
				}while(fabs(d)>EPS && n<NMAX);

				if(n==NMAX){
					fprintf(fp,"%f fault\n",temperature);
				}else{
					fprintf(fp,"%f %f\n",temperature,fabs(mag));
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

double RightHandSide(double mag,double temperature,double c){
 
    int k,l;
    double summation_k;
    double summation_l;
 
    summation_l = 0.0;
 
    for(l=0;l<=30;l++){
 
        summation_k = 0.0;
 
        for(k=0;k<=l;k++){
			if(temperature>0){
				summation_k += ( factorial(l)/( factorial(k)*factorial(l-k) ) )
								*( pow(pow((1+mag)/2.0,2)+pow((1-mag)/2.0,2),k)*pow(0.5*(1+mag)*(1-mag),l-k)/( 1.0+exp(-2.0*(2.0*k-l)/temperature) ) );
			}else if(temperature==0.0){
				if(2*k-l>0){
					summation_k += ( factorial(l)/( factorial(k)*factorial(l-k) ) )
					*( pow(pow((1+mag)/2.0,2)+pow((1-mag)/2.0,2),k)*pow(0.5*(1+mag)*(1-mag),l-k) );
				}else if(2*k-l==0){
					summation_k += 0.5*( factorial(l)/( factorial(k)*factorial(l-k) ) )
					*( pow(pow((1+mag)/2.0,2)+pow((1-mag)/2.0,2),k)*pow(0.5*(1+mag)*(1-mag),l-k) );
				}else if(2*k-l<0){
					summation_k += 0.0;
				}
			}

		}
 
        summation_l += ( exp(-c)*pow(c,l)/factorial(l) )*summation_k;
    }
 
    return -(mag+1.0)+2.0*summation_l;
 
}

double dRightHandSide(double mag,double temperature,double c){
	
	int k,l;
	double summation_k,summation_l;

	summation_l = 0.0;
	
	for(l=0;l<=30;l++){

		summation_k = 0.0;

		for(k=0;k<=l;k++){

			if(temperature>0){
				summation_k += ( factorial(l)/(factorial(k)*factorial(l-k)) )*
					( pow( pow((1+mag)/2.0,2)+pow((1-mag)/2.0,2),k-1 ) )*( pow((1-pow(mag,2))/2.0,l-k-1) )*( k-l*(1+pow(mag,2))/2.0 )
					/( 1+exp(-2*(2*k-l)/temperature) );
			}else if(temperature==0.0){
				if(2*k-l>0){
					summation_k += ( factorial(l)/(factorial(k)*factorial(l-k)) )*
									( pow( pow((1+mag)/2.0,2)+pow((1-mag)/2.0,2),k-1 ) )*( pow((1-pow(mag,2))/2.0,l-k-1) )*( k-l*(1+pow(mag,2))/2.0 );
				}else if(2*k-l==0){
					summation_k += 0.5*( factorial(l)/(factorial(k)*factorial(l-k)) )*
									( pow( pow((1+mag)/2.0,2)+pow((1-mag)/2.0,2),k-1 ) )*( pow((1-pow(mag,2))/2.0,l-k-1) )*( k-l*(1+pow(mag,2))/2.0 );
				}else if(2*k-l<0){
					summation_k += 0.0;
				}
			}

		summation_l += ( exp(-c)*pow(c,l)/factorial(l) )*summation_k;

		}
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