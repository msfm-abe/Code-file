#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include "dSFMT.h"

#define SIZE 16384
#define EDGE SIZE*c
#define size (int)((SIZE*c)*0.5)
#define MCS 100000
#define SAMPLE 80//sample数（グラフの数）
#define TEMPERATURE 0.5
#define IM 1.0

dsfmt_t dsfmt;

int main(void){

	int i,j,k;
	int seed;
	int *spin;
	double c;
	double magnetization;
	int **List;
	int *deg;
	int choice;
	int step;
	int sample;
	double probability;
	int hi;
	double mcurrent;
	double *memory;
	double ave_magnetization;
	double ave_magnetization2;
	double denominator;

	FILE *fp;
	fp = fopen("2-body_c-Magnetization_transition_c=0.5-3.5_T=0.5_m0=1.0.dat","wb");

	spin = (int *)malloc(sizeof(int)*SIZE);
	deg = (int *)malloc(sizeof(int)*SIZE);
	memory = (double *)malloc(sizeof(double)*SAMPLE);


	for(c=0.5;c<=3.5;c+=0.1){

		for(i=0;i<SAMPLE;i++){
			memory[i] = 0.0;
		}

		List = (int **)malloc(sizeof(int *)*SIZE);
		for(i=0;i<SIZE;i++){
			List[i] = (int *)malloc(sizeof(int)*size);
		}

		for(sample=0;sample<SAMPLE;sample++){

			seed = (unsigned int)time(NULL);
			dsfmt_init_gen_rand(&dsfmt,seed);

			for(i=0;i<SIZE;i++){
				for(j=0;j<size;j++){
					List[i][j] = -1;
				}
				deg[i] = 0;
			}

			//initial value for spins
			for(i=0;i<SIZE;i++){
				spin[i] = ( (double)dsfmt_genrand_close_open(&dsfmt)<=(IM+1.0)/2.0 )? 1:-1;
			}

			//create the graph

			for(i=0;i<EDGE;i++){
				choice = (int)((double)dsfmt_genrand_close_open(&dsfmt)*SIZE);

				if(i!=choice){
					List[choice][deg[choice]] = i;
					deg[choice]++;
				}else{
					i--;
					continue;
				}
			}

			magnetization = 0.0;

			for(i=1;i<=MCS;i++){

				for(step=0;step<SIZE;step++){
					choice = (int)(dsfmt_genrand_close_open(&dsfmt)*(SIZE));

					hi = 0;
					for(j=0;j<deg[choice];j+=1){
						hi += spin[List[choice][j]];
					}

					if(TEMPERATURE>0.0){
						probability = 1.0/( 1.0+exp(-2.0*hi/TEMPERATURE) );
					}else if(TEMPERATURE==0.0){
						if(hi>0){
							probability = 1.0;
						}else if(hi==0){
							probability = 0.5;
						}else if(hi<0){
							probability = 0.0;
						}
					}

					if( (double)dsfmt_genrand_close_open(&dsfmt)<=probability ){
						spin[choice] = 1;
					}else{
						spin[choice] = -1;
					}
				}//step

				if(i>1000){
					mcurrent = 0.0;
					for(k=0;k<SIZE;k++){
						mcurrent += (double)spin[k];
					}
					mcurrent /= SIZE;

					magnetization += mcurrent;
				}
		
			}//mcs

			denominator = (double)(MCS-1000);
			memory[sample] = magnetization/denominator;

		}//sample

		ave_magnetization = 0.0;
		ave_magnetization2 = 0.0;

		for(i=0;i<SAMPLE;i++){
			ave_magnetization += memory[i];
			ave_magnetization2 += pow(memory[i],2);
		}
		ave_magnetization /= SAMPLE;
		ave_magnetization2 /= SAMPLE;
	
		fprintf(fp,"%f %f %f\n",c,fabs(ave_magnetization),sqrt( ave_magnetization2-pow(ave_magnetization,2.0) )/sqrt( (double)SAMPLE ));


		for(i=0;i<SIZE;i++){
			free((void *)List[i]);
		}
		free((void *)List);
		
	}//c

	free((void *)spin);
	free((void *)deg);
	free((void *)memory);

	fclose(fp);

	//getchar();
	return 0;
}