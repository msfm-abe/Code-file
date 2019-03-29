#include<stdlib.h>

int *ikakuho1(int size){

	int *p;
	int i;

	p = (int *)malloc(sizeof(int)*size);
	if(p==NULL) exit(0);

	for(i=0;i<size;i++){
		p[i] = 0;
	}

	return p;

}

double *dkakuho1(int size){

	double *p;
	int i;

	p = (double *)malloc(sizeof(int)*size);
	if(p==NULL) exit(0);

	for(i=0;i<size;i++){
		p[i] = 0.0;
	}

	return p;

}

int **ikakuho2(int row,int column){

	int **p;
	int i,j;

	p = (int **)malloc(sizeof(int *)*row);
	if(p==NULL) exit(0);
	for(i=0;i<row;i++){
		p[i] = (int *)malloc(sizeof(int)*column);
		if(p[i]==NULL) exit(0);
	}

	for(i=0;i<row;i++){
		for(j=0;j<column;j++){
			p[i][j] = 0;
		}
	}

	return p;

}

double **dkakuho2(int row,int column){

	double **p;
	int i,j;

	p = (double **)malloc(sizeof(double *)*row);
	if(p==NULL) exit(0);
	for(i=0;i<row;i++){
		p[i] = (double *)malloc(sizeof(double)*column);
		if(p[i]==NULL) exit(0);
	}

	for(i=0;i<row;i++){
		for(j=0;j<column;j++){
			p[i][j] = 0.0;
		}
	}

	return p;

}

void ikaiho1(int *f){

	free(f);

}

void dkaiho1(double *f){
	
	free(f);

}

void ikaiho2(int **f,int row){

	int i;
	for(i=0;i<row;i++){
		free(f[i]);
	}
	free(f);

}

void dkaiho2(double **f,int row){

	int i;
	for(i=0;i<row;i++){
		free(f[i]);
	}
	free(f);

}