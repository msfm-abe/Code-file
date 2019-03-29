#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"dSFMT.h"

#define SIZE 512
#define size (int)SIZE*0.5
#define SAMPLE 800
#define MCS 100001
#define IM 0.0
#define EDGE (int)SIZE*c
//#define TEMPERATURE 0.5
//#define C 3.74

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
	double c;//=C;
	double temperature;//=TEMPERATURE;
	double probability;
	double **hairetsu;
	double correlation=0.0;
	double correlation2=0.0;
	char filename[100];
	int ketasuu,kakunin;
	double Sample=(double)SAMPLE;

	FILE *fp;

	seed = (unsigned int)time(NULL);
	dsfmt_init_gen_rand(&dsfmt,seed);


	spin_0 = (int *)malloc(sizeof(int)*SIZE);

	spin = (int *)malloc(sizeof(int)*SIZE);

	List = (int **)malloc(sizeof(int *)*SIZE);
	for(i=0;i<SIZE;i++){
		List[i] = (int *)malloc(sizeof(int)*SIZE);
	}

	deg = (int *)malloc(sizeof(int)*SIZE);

	hairetsu = (double **)malloc(sizeof(double *)*MCS);
	for(i=0;i<MCS;i++){
		hairetsu[i] = (double *)malloc(sizeof(double)*SIZE);
	}

	
	for(c=4.2;c<=4.6;c+=0.2){

		//sprintf(filename,"3-body_UndirectedVer_glass_c-T_c=%3.1f.dat",c);
		//fp = fopen(filename,"wb");

		//fprintf(fp,"#c,T,MCS,correlation,errorbar\n");

		pair = (int **)malloc(sizeof(int *)*EDGE);
		for(i=0;i<EDGE;i++){
			pair[i] = (int *)malloc(sizeof(int)*2);
		}

		for(temperature=0.0;temperature<=1.0;temperature+=0.1){

			sprintf(filename,"3-body_UndirectedVer_glass_c-T_c=%3.1f_T=%3.1f.dat",c,temperature);
			fp = fopen(filename,"wb");

			for(i=0;i<MCS;i++){
				for(j=0;j<SIZE;j++){
					hairetsu[i][j] = 0.0;
				}
			}

			for(sample=0;sample<SAMPLE;sample++){

				for(i=0;i<SIZE;i++){

					spin[i] = ( (double)dsfmt_genrand_close_open(&dsfmt)<=(IM+1)/2.0 )? 1:-1;

				}

				for(i=0;i<SIZE;i++){
					for(j=0;j<SIZE;j++){
						List[i][j] = -1;
					}
					deg[i] = 0;
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

						List[pair[i][0]][deg[pair[i][0]]] = choice;
						deg[pair[i][0]]++;
						List[pair[i][0]][deg[pair[i][0]]] = pair[i][1];
						deg[pair[i][0]]++;

						List[pair[i][1]][deg[pair[i][1]]] = choice;
						deg[pair[i][1]]++;
						List[pair[i][1]][deg[pair[i][1]]] = pair[i][0];
						deg[pair[i][1]]++;

					}else{

						--i;
						continue;

					}

				}

				for(mcs=1;mcs<MCS;mcs++){

					for(step=0;step<SIZE;step++){

						hi = 0;
						choice = (int)(dsfmt_genrand_close_open(&dsfmt)*(SIZE));

						for(i=0;i<deg[choice];i+=2){
							hi += spin[choice]*spin[List[choice][i]]*spin[List[choice][i+1]];
						}

						if(temperature>0.0){
							probability = 1.0/( 1.0+exp(2*hi/temperature) );
						}else if(temperature==0.0){
							if(hi>0){
								probability = 0.0;
							}else if(hi==0){
								probability = 0.5;
							}else if(hi<0){
								probability = 1.0;
							}
						}

						spin[choice] = ( (double)dsfmt_genrand_close_open(&dsfmt)<=probability )? 1:-1;

					}//step

					if(mcs==1000){

						for(i=0;i<SIZE;i++){
							spin_0[i] = spin[i];
						}

					}

					if(mcs>=1000){

						for(i=0;i<SIZE;i++){

							hairetsu[mcs][i] += (double)spin_0[i]*spin[i];

						}

					}

				}//mcs


			}//sample


			for(mcs=0;mcs<MCS;mcs++){
				for(i=0;i<SIZE;i++){

					hairetsu[mcs][i] /= Sample;

				}
			}


			for(mcs=0;mcs<MCS;mcs++){

				if(mcs>=1000){

					ketasuu = 0;
					kakunin = mcs;
				
					while(kakunin!=0){
						kakunin /= 10;
						++ketasuu;
					}

					if( (mcs%(int)pow(10.0,ketasuu-1))==0 ){

						correlation = 0.0;
						correlation2 = 0.0;

						for(i=0;i<SIZE;i++){
							correlation += hairetsu[mcs][i];
							correlation2 += pow(hairetsu[mcs][i],2);
						}

						correlation /= SIZE;
						correlation2 /= SIZE;

						fprintf(fp,"%f %f %6d %f %f\n",c,temperature,mcs,correlation,sqrt( correlation2-pow(correlation,2) )/sqrt( Sample ));

					}

				}

			}//second_mcs

			//fprintf(fp,"\n");
			fclose(fp);

		}//temperature


		for(i=0;i<EDGE;i++){
			free((void *)pair[i]);
		}
		free((void *)pair);

		//fclose(fp);

	}//c


	free((void *)spin_0);
	free((void *)spin);
	free((void *)deg);
	for(i=0;i<SIZE;i++){
		free((void *)List[i]);
	}
	free((void *)List);
	for(i=0;i<MCS;i++){
		free((void *)hairetsu[i]);
	}
	free((void *)hairetsu);

	//getchar();
	return 0;
}
