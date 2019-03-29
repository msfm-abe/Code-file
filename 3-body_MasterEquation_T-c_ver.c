#include<stdio.h>
#include<stdlib.h>
#include<math.h>
 
double factorial(int k);
double RightHandSide(double temp,double c);
 
#define mag 1.0
 
int main(){
 
    double temp;
    double c;
    char filename[100];
 
    FILE *fp;
	fp = fopen("testfile20160908.dat","wb");

	for(c=3.0;c<=5.0;c+=0.001){
		sprintf(filename,"3-body_MasterEquation_T-c_ver_m0=1.0_c=%f.dat",c);

 
    return 0;
 
}
 
double RightHandSide(double temp,double c){
 
    int k,l;
    double summation_k;
    double summation_l;
 
    summation_l = 0.0;
 
    for(l=0;l<=20;l++){
 
        summation_k = 0.0;
 
        for(k=0;k<=l;k++){
            summation_k += ( factorial(l)/( factorial(k)*factorial(l-k) ) )
                *( pow(pow((1+mag)/2.0,2)+pow((1-mag)/2.0,2),k)*pow(0.5*(1+mag)*(1-mag),l-k)/( 1.0+exp(-2.0*(2.0*k-l)/temp) ) );
        }
 
        summation_l += ( exp(-c)*pow(c,l)/factorial(l) )*summation_k;
    }
 
    return -(mag+1.0)+2.0*summation_l;
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