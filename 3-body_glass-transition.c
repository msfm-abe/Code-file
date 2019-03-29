#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"dSFMT.h"

#define SIZE 16384
#define size (int)SIZE*0.5
#define SAMPLE 80
#define MCS 100000
#define IM 0.0
#define EDGE (int)SIZE*c
#define TEMPERATURE 0.5
#define C 3.0

dsfmt_t dsfmt;

int main(void){

	int i,j,k;
	int seed;
	int choice;
	int mcs;
	int step;
	int hi;
	int sample;
	int *spin_0;
	int *spin;
	int *deg;
	int **List;
	int **pair;
	double c=C;
	double temperature=TEMPERATURE;
	double probability;
	double **hairetsu;
	double *heikin;
	double correlation=0.0;

	FILE *fp;
	fp = fopen("glasstest.dat","wb");

	seed = (unsigned int)time(NULL);
    dsfmt_init_gen_rand(&dsfmt,seed);

	
	spin_0 = (int *)malloc(sizeof(int)*SIZE);

	spin = (int *)malloc(sizeof(int)*SIZE);

	pair = (int **)malloc(sizeof(int *)*EDGE);
	for(i=0;i<EDGE;i++){
		pair[i] = (int *)malloc(sizeof(int)*2);
	}

	List = (int **)malloc(sizeof(int *)*SIZE);
	for(i=0;i<SIZE;i++){
		List[i] = (int *)malloc(sizeof(int)*size);
	}

	deg = (int *)malloc(sizeof(int)*SIZE);

	hairetsu = (double **)malloc(sizeof(double *)*SAMPLE);
	for(i=0;i<SAMPLE;i++){
		hairetsu[i] = (double *)malloc(sizeof(double)*SIZE);
	}

	heikin = (double *)malloc(sizeof(double)*SIZE);


	for(sample=0;sample<SAMPLE;sample++){

		for(i=0;i<SIZE;i++){

			spin[i] = ( (double)dsfmt_genrand_close_open(&dsfmt)<=(IM+1)/2.0 )? 1:-1;

		}

		for(i=0;i<EDGE;i++){
			for(j=0;j<2;j++){
				pair[i][j] = -1;
			}
		}

		for(i=0;i<SIZE;i++){
			for(j=0;j<size;j++){
				List[i][j] = -1;
			}
			deg[i] = 0;
		}

		for(i=0;i<SAMPLE;i++){
			for(j=0;j<SIZE;j++){
				hairetsu[i][j] = 0.0;
			}
		}

		for(i=0;i<SIZE;i++){
			
			heikin[i] = 0.0;

		}


		for(i=0;i<EDGE;i++){

				pair[i][0] = (int)(dsfmt_genrand_close_open(&dsfmt)*(SIZE));
				pair[i][1] = (int)(dsfmt_genrand_close_open(&dsfmt)*(SIZE));

				if(pair[i][0]==pair[i][1]){
					i--;
					continue;
				}
		}

		for(i=0;i<EDGE;i++){

				choice = (int)(dsfmt_genrand_close_open(&dsfmt)*(SIZE));

				if(choice!=pair[i][0] && choice!=pair[i][1]){
					List[choice][deg[choice]] = pair[i][0];
					deg[choice]++;
					List[choice][deg[choice]] = pair[i][1];
					deg[choice]++;
				}else{
					i--;
					continue;
				}
		}


		for(mcs=1;mcs<=MCS;mcs++){

			if(mcs==1000){

				for(i=0;i<SIZE;i++){
					spin_0[i] = spin[i];
				}

			}

			for(step=0;step<SIZE;step++){

				hi = 0;
				choice = (int)(dsfmt_genrand_close_open(&dsfmt)*(SIZE));

				for(i=0;i<deg[choice];i+=2){
					hi += spin[List[choice][i]]*spin[List[choice][i+1]];
				}

				if(temperature>0.0){
					probability = 1.0/( 1.0+exp(-2*hi/TEMPERATURE) );
				}else if(temperature==0.0){
					if(hi>0){
						probability = 1.0;
					}else if(hi==0){
						probability = 0.5;
					}else if(hi<0){
						probability = 0.0;
					}
				}

				spin[choice] = ( (double)dsfmt_genrand_close_open(&dsfmt)<=probability )? 1:-1;

			}//step

		}//mcs

		for(i=0;i<SIZE;i++){

			hairetsu[sample][i] = spin_0[i]*spin[i];

		}

	}//sample

	for(i=0;i<SIZE;i++){
		for(j=0;j<SAMPLE;j++){
			heikin[i] += hairetsu[j][i];
		}
		heikin[i] /= SAMPLE;//‚±‚ÌŽž“_‚Åheikin‚ÌŠe¬•ª‚É‚Í<s(0)*s(\infty)>‚ª“ü‚Á‚Ä‚¢‚é
	}

	for(i=0;i<SIZE;i++){
		correlation += heikin[i];
	}

	correlation /= SIZE;

	fprintf(fp,"%f %f %f\n",c,temperature,correlation);

	fclose(fp);


	free((void *)spin_0);
	free((void *)spin);
	free((void *)deg);
	for(i=0;i<SIZE;i++){
		free((void *)List[i]);
	}
	free((void *)List);
	for(i=0;i<EDGE;i++){
		free((void *)pair[i]);
	}
	free((void *)pair);

	//getchar();
	return 0;
}