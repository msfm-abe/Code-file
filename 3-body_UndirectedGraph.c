#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"dSFMT.h"

#define SIZE 20
#define EDGE (int)c*SIZE
#define C 1.0

dsfmt_t dsfmt;

int main(void){

	int seed;
	int i,j;
	int **List;
	int *deg;
	int **pair;
	double c=C;
	int choice;
	FILE *fp;
	
	fp = fopen("mukougurahu.dat","wb");

	seed = (unsigned int)time(NULL);
	dsfmt_init_gen_rand(&dsfmt,seed);

	List = (int **)malloc(sizeof(int *)*SIZE);
	for(i=0;i<SIZE;i++){
		List[i] = (int *)malloc(sizeof(int)*SIZE);
	}

	pair = (int **)malloc(sizeof(int *)*EDGE);
	for(i=0;i<EDGE;i++){
		pair[i] = (int *)malloc(sizeof(int)*2);
	}

	deg = (int *)malloc(sizeof(int)*SIZE);

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

	for(i=0;i<SIZE;i++){

		fprintf(fp,"deg[%2d]=%2d ",i,deg[i]);

		for(j=0;j<SIZE;j++){

			fprintf(fp,"%3d ",List[i][j]);

		}
		fprintf(fp,"\n");
	}

	fclose(fp);

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