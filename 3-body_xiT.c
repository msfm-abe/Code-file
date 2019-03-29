#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"dSFMT.h"
 
#define SIZE 10
#define C 3.8
#define EDGE (int)(SIZE*c)
#define size (int)((SIZE*c)*0.5)
#define MCS 5000
#define SAMPLE 2
#define GRAPH 3
#define InitialMag 1.0
#define TEMPERATURE 0.5
 
dsfmt_t dsfmt;
 
int main(void){

	int i,j,k;
    int seed;
    int *spin;
    double c=C;
    double magnetization;
    int **List;
    int **pair;
    int *deg;
    int choice;
    int step;
    int mcs;
    int sample;
    int graph;
    double probability;
    int hi;
    double mcurrent;
    double **memory;
	double tr_average1;
	double tr_average2;
    double **tr_average;
    double xiT;
	double xiT2;

	FILE *fp;
    fp = fopen("xiT_test.dat","wb");
     
    seed = (unsigned int)time(NULL);
    dsfmt_init_gen_rand(&dsfmt,seed);

	spin = (int *)malloc(sizeof(int)*SIZE);

	memory = (double **)malloc(sizeof(double *)*SAMPLE);
	for(i=0;i<SAMPLE;i++){
		memory[i] = (double *)malloc(sizeof(double)*MCS);
	}

	tr_average = (double **)malloc(sizeof(double *)*GRAPH);
	for(i=0;i<GRAPH;i++){
		tr_average[i] = (double *)malloc(sizeof(double)*MCS);
	}

	for(i=0;i<SAMPLE;i++){
		for(j=0;j<MCS;j++){
			memory[i][j] = 0.0;
		}
	}

	//for(i=1;i<=7;i++){
	//	printf("exp(-4.0*%d)=%f\n",i,exp(-4.0*i));
	//}

	//for(i=0;i<SAMPLE;i++){
	//	for(j=0;j<MCS;j++){
	//		printf("memory[%d][%d]=%f\n",i,j,memory[i][j]);
	//	}
	//}

	for(i=0;i<GRAPH;i++){
		for(j=0;j<MCS;j++){
			tr_average[i][j] = 0.0;
		}
	}

	//for(i=0;i<GRAPH;i++){
	//	for(j=0;j<MCS;j++){
	//		printf("tr_average[%d][%d]=%f\n",i,j,tr_average[i][j]);
	//	}
	//}


	for(graph=0;graph<GRAPH;graph++){

		pair = (int **)malloc(sizeof(int *)*EDGE);
		for(i=0;i<EDGE;i++){
			pair[i] = (int *)malloc(sizeof(int)*2);
		}

		deg = (int *)malloc(sizeof(int)*SIZE);

		List = (int **)malloc(sizeof(int *)*SIZE);
		for(i=0;i<SIZE;i++){
			List[i] = (int *)malloc(sizeof(int)*size);
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

		for(i=0;i<EDGE;i++){

			pair[i][0] = (int)(dsfmt_genrand_close_open(&dsfmt)*(SIZE));
			pair[i][1] = (int)(dsfmt_genrand_close_open(&dsfmt)*(SIZE));

			if(pair[i][0]==pair[i][1]){
				i--;
				continue;
			}
		}

		//for(i=0;i<EDGE;i++){
		//	printf("pair[%d]  %d %d\n",i,pair[i][0],pair[i][1]);
		//}

		for(i=0;i<EDGE;i++){

			choice =  (int)(dsfmt_genrand_close_open(&dsfmt)*(SIZE));

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

		//for(i=0;i<SIZE;i++){
		//	printf("L[%d] deg[%d]=%d ",i,i,deg[i]);
		//	for(j=0;j<size;j++){
		//		printf("%d ",List[i][j]);
		//	}
		//	printf("\n");
		//}

		for(sample=0;sample<SAMPLE;sample++){

			//initialization for spins
			for(i=0;i<SIZE;i++){
				if( (double)dsfmt_genrand_close_open(&dsfmt)<=(InitialMag+1.0)/2.0 ){
					spin[i] = 1;
				}else{
					spin[i] = -1;
				}
			}

			magnetization = 0.0;

			for(mcs=1;mcs<=MCS;mcs++){

				for(step=0;step<SIZE;step++){
					choice = (int)(dsfmt_genrand_close_open(&dsfmt)*(SIZE));

					hi = 0;
					for(i=0;i<deg[choice];i+=2){
						hi += spin[List[choice][i]]*spin[List[choice][i+1]];
					}
					//printf("hi=%d\n",hi);

					probability = 1.0/( 1.0+exp(-2.0*hi/TEMPERATURE) );
					//printf("%f %f\n",probability,exp(-2.0*hi/TEMPERATURE));

					if( (double)dsfmt_genrand_close_open(&dsfmt)<=probability ){
						spin[choice] = 1;
					}else{
						spin[choice] = -1;
					}

				}//step

				mcurrent = 0;
				for(i=0;i<SIZE;i++){
					mcurrent += (double)spin[i];
				}
				mcurrent /= (double)SIZE;

				magnetization += mcurrent;

				memory[sample][mcs-1] = magnetization/mcs;

			}//mcs

			//for(i=0;i<MCS;i++){
			//	printf("m[%d][%d]=%f ",sample,i,memory[sample][i]);
			//}
			//printf("\n");

		}//sample

		//calculate <m> and <m^2>

		
		for(i=0;i<MCS;i++){
			tr_average1 = 0.0;
			tr_average2 = 0.0;
			for(j=0;j<SAMPLE;j++){
				tr_average1 += memory[j][i];
				tr_average2 += pow(memory[j][i],2);
			}
			tr_average1 /= SAMPLE;//<m>
			tr_average2 /= SAMPLE;//<m^2>

			tr_average[graph][i] = tr_average2-pow(tr_average1,2);
		}

	}//graph

	//for(i=0;i<GRAPH;i++){
	//	for(j=0;j<MCS;j++){
	//		printf("t_a[%d][%d]=%f ",i,j,tr_average[i][j]);
	//	}
	//	printf("\n\n");
	//}

	for(mcs=0;mcs<MCS;mcs++){
		xiT = 0.0;
		xiT2 = 0.0;
		//if( (mcs%=1000)==999 ){
			for(graph=0;graph<GRAPH;graph++){
				xiT += tr_average[graph][mcs];
				xiT2 += pow(tr_average[graph][mcs],2);
			}
			xiT /= GRAPH;
			xiT2 /= GRAPH;

			//xiT *= SIZE;

			fprintf(fp,"%d %f %f\n",mcs+1,fabs(SIZE*xiT),sqrt( xiT2-pow(xiT,2) )/sqrt( (double)GRAPH ));
		//}
	}

	fclose(fp);

	free((void *)spin);
	for(i=0;i<SAMPLE;i++){
		free((void *)memory[i]);
	}
	free((void *)memory);
	for(i=0;i<GRAPH;i++){
		free((void *)tr_average[i]);
	}
	free((void *)tr_average);
	for(i=0;i<EDGE;i++){
		free((void *)pair[i]);
	}
	free((void *)pair);
	for(i=0;i<SIZE;i++){
		free((void *)List[i]);
	}
	free((void *)List);
	free((void *)deg);

	return 0;

}