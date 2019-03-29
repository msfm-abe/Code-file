#include<cstdlib>

int *intArrya1(int row){	
	int *p = 0;
	p = new int[row];
	if(p==0) exit(0);

	for(int i=0;i<row;i++)
		p[i] = 0;

	return p;
}

int **intArray2(int row,int column){
	int **p = 0;
	p = new int*[row];
	if(p==0) exit(0);
	for(int i=0;i<row;i++){
		p[i] = new int[column];
		if(p[i]==0) exit(0);
	}

	for(int i=0;i<row;i++){
		for(int j=0;j<column;j++){
			p[i][j] = 0;
		}
	}

	return p;
}

double *doubleArray1(int row){
	double *p = 0;
	p = new double[row];
	if(p==0) exit(0);
	
	for(int i=0;i<row;i++)
		p[i] = 0.0; 

	return p;
}

double **doubleArray2(int row,int column){
	double **p = 0;
	p = new double*[row];
	if(p==0) exit(0);
	for(int i=0;i<row;i++){
		p[i] = new double[column];
		if(p[i]==0) exit(0);
	}

	for(int i=0;i<row;i++){
		for(int j=0;j<column;j++){
			p[i][j] = 0.0;
		}
	}

	return p;
}

void ikaiho1(int *p){
	delete p;
}

void ikaiho2(int **p,int row){
	for(int i=0;i<row;i++)
		delete p[i];

	delete p;
}

void dkaiho1(double *p){
	delete p;
}

void dkaiho2(double **p,int row){
	for(int i=0;i<row;i++)
		delete p[i];

	delete p;
}