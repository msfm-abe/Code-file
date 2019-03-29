#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include "dSFMT.h"

#define SIZE 512
#define EDGE SIZE*c
#define size (int)((SIZE*c)*0.5)
#define MCS 100000
#define SAMPLE 800//sample数（グラフの数）
#define temperature 0.5

dsfmt_t dsfmt;

int main(void){

	int i,j,k;
	int seed;
	int *spin;
	double c;
	double magnetization;
	int **List;
	int **pair;
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
	fp = fopen("3-body_test2.dat","wb");

	spin = (int *)malloc(sizeof(int)*SIZE);
	deg = (int *)malloc(sizeof(int)*SIZE);
	memory = (double *)malloc(sizeof(double)*SAMPLE);

	for(i=0;i<SAMPLE;i++){
		memory[i] = 0.0;
	}


	for(c=3.0;c<=5.0;c+=0.1){

		List = (int **)malloc(sizeof(int *)*SIZE);
		for(i=0;i<SIZE;i++){
			List[i] = (int *)malloc(sizeof(int)*size);
		}
	
		pair = (int **)malloc(sizeof(int *)*EDGE);
		for(i=0;i<EDGE;i++){
			pair[i] = (int *)malloc(sizeof(int)*2);
		}


		for(sample=0;sample<SAMPLE;sample++){

			for(i=0;i<SIZE;i++){
				for(j=0;j<size;j++){
					List[i][j] = -1;
				}
				deg[i] = 0;
			}

			for(i=0;i<EDGE;i++){
				for(j=0;j<2;j++){
					pair[i][j] = -1;
				}
			}

			//initial value for spins
			for(i=0;i<SIZE;i++){
				spin[i] = 1;
			}

			//create the graph
	
			seed = (unsigned int)time(NULL);
			dsfmt_init_gen_rand(&dsfmt,seed);

			for(i=0;i<EDGE;i++){
				pair[i][0] = (int)(dsfmt_genrand_close_open(&dsfmt)*(SIZE));
				pair[i][1] = (int)(dsfmt_genrand_close_open(&dsfmt)*(SIZE));

				if(pair[i][0]==pair[i][1]){
					i--;
					continue;
				}
			}

			for(i=0;i<EDGE;i++){
				choice = (int)((double)dsfmt_genrand_close_open(&dsfmt)*SIZE);
		
				if(choice!=pair[i][0] && choice!=pair[i][1]){
					List[choice][deg[choice]] = pair[i][0];
					deg[choice] += 1;
					List[choice][deg[choice]] = pair[i][1];
					deg[choice] += 1;
 				}else{
					i--;
					continue;
				}
			}

			//for(i=0;i<EDGE;i++){
			//	printf("pair[%d][0]=%d pair[%d][1]=%d\n",i,pair[i][0],i,pair[i][1]);
			//}
			//printf("\n");
			//for(i=0;i<SIZE;i++){
			//	printf("deg[%d]=%d  ",i,deg[i]);
			//	for(j=0;j<size;j++){
			//		printf("%d ",List[i][j]);
			//	}
			//	printf("\n");
			//}

			magnetization = 0.0;

			for(i=1;i<=MCS;i++){

				for(step=0;step<SIZE;step++){
					choice = (int)(dsfmt_genrand_close_open(&dsfmt)*(SIZE));

					hi = 0;
					for(j=0;j<deg[choice];j+=2){
						hi += spin[List[choice][j]]*spin[List[choice][j+1]];
					}
					probability = 1.0/( 1.0+exp(-2.0*hi/temperature) );

					//printf("hi=%d probability=%f\n",hi,probability);

					if( (double)dsfmt_genrand_close_open(&dsfmt)<probability ){
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
		for(i=0;i<EDGE;i++){
			free((void *)pair[i]);
		}
		free((void *)pair);
		
	}//c

	//for(i=0;i<SAMPLE;i++){
	//	printf("memory[%d]=%f\n",i,memory[i]);
	//}

	free((void *)spin);
	free((void *)deg);
	free((void *)memory);

	fclose(fp);

	//getchar();
	return 0;
}