#include<stdio.h>
#include<stdlib.h>
#include<math.h>
 
#define temperature 0.5
#define InitialMagnetization 1.0
 
double RightHandSide(double mag,double c);
double factorial(int k);
 
int main(void){
 
    int i,j,k;
    double mag;
    double k1,k2,k3,k4;
    double h,a,b,n;
    double time;
    double c;
	double im;
	char filename[100];
 
    FILE *fp;

   //fp = fopen("3-body_MasterEquation_M.ver_c_izon_T=0.5_InitialMagnetization=1.0.dat","wb");
 
    a = 0.0;
    b = 1000.0;
    n = 10000.0;
    h = (b-a)/n;//h is fixed to 0.1;

	im = 0.6;
 
	//for(im=0.8;im<=InitialMagnetization;im+=0.1){
		sprintf(filename,"3-body_MasterEquation_M.ver_c_izon_T=0.5_InitialMagnetization=%3.1f.dat",im);
		fp = fopen(filename,"wb");
 
		for(c=3.0;c<=5.0;c+=0.01){

			mag = im;
 
			for(time=a;time<=b;time+=h){
				k1 = RightHandSide(mag,c);
				k2 = RightHandSide(mag+h*k1/2.0,c);
				k3 = RightHandSide(mag+h*k2/2.0,c);
				k4 = RightHandSide(mag+h*k3,c);
 
				mag = mag + h*(k1+2.0*k2+2.0*k3+k4)/6.0;
			}
 
			fprintf(fp,"%f %f\n",c,fabs(mag));
		}
 
		fclose(fp);

	//}
 
    return 0;
}
 
double RightHandSide(double mag,double c){
 
    int k,l;
    double summation_k;
    double summation_l;
 
    summation_l = 0.0;
 
    for(l=0;l<=20;l++){
 
        summation_k = 0.0;
 
        for(k=0;k<=l;k++){
            summation_k += ( factorial(l)/( factorial(k)*factorial(l-k) ) )
                *( pow(pow((1+mag)/2.0,2)+pow((1-mag)/2.0,2),k)*pow(0.5*(1+mag)*(1-mag),l-k)/( 1.0+exp(-2.0*(2.0*k-l)/temperature) ) );

			//if( (2.0*k-l)>0 ){
			//	summation_k += ( factorial(l)/( factorial(k)*factorial(l-k) ) )
   //             *( pow(pow((1+mag)/2.0,2)+pow((1-mag)/2.0,2),k)*pow(0.5*(1+mag)*(1-mag),l-k) );
			//}else if( (2.0*k-l)==0 ){
			//	summation_k += 0.5*( factorial(l)/( factorial(k)*factorial(l-k) ) )
   //             *( pow(pow((1+mag)/2.0,2)+pow((1-mag)/2.0,2),k)*pow(0.5*(1+mag)*(1-mag),l-k) );
			//}else if( (2.0*k-l)<0 ){
			//	summation_k += 0.0;
			//}
		}
 
        summation_l += ( pow(c,l)/factorial(l) )*summation_k;
    }
 
    return -(mag+1.0)+2.0*exp(-c)*summation_l;
 
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