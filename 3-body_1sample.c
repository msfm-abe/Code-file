#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include "dSFMT.h"

#define SIZE 16384
#define size (int)((SIZE*c)*0.2)
#define MCS 100001
#define SAMPLE 1//sample数（グラフの数）

int STEP;

void get_average(double **);

dsfmt_t dsfmt;

int main(void)
{
    int mcs;
	int sample;
    int i;
	int j;
    int i1;
	int cite;
	int seed;
    double temperature;
    double denominator;
	double mcurrent;
	double magnetization;
	double probability;
    double r;
    double c;
	int **List;
	double **Sample;
	int *deg;
	int *spin;
	int *hi;
	int **pair;

    c = 3.9;
	temperature = 0.5;

	seed = (unsigned int)time(NULL);
	dsfmt_init_gen_rand(&dsfmt, seed);

	Sample = (double **)malloc(sizeof(double *)*MCS);

	for(i=0;i<MCS;i++)
	{
		Sample[i] = (double *)malloc(sizeof(double)*SAMPLE);
	}

	pair = (int **)malloc(sizeof(int *)*((int)SIZE*c));

	for(i=0;i<((int)SIZE*c);i++)
	{
		pair[i] = (int *)malloc(sizeof(int)*2);
	}

	List = (int **)malloc(sizeof(int *)*SIZE);
	
	for(i=0;i<SIZE;i++)
	{
		List[i] = (int *)malloc(sizeof(int)*size);
	}

	deg = (int *)malloc(sizeof(int)*SIZE);
	
	hi = (int *)malloc(sizeof(int)*SIZE);
	
	spin = (int *)malloc(sizeof(int)*SIZE);

	for(STEP=1;STEP<=2;STEP++){

		for(i=0;i<MCS;i++)
		{
			for(j=0;j<SAMPLE;j++)
			{
				Sample[i][j] = 0.0;
			}
		}
	
		for(sample=0;sample<SAMPLE;sample++)
		{

			for(i=0;i<SIZE;i++)
			{
				for(j=0;j<size;j++)
				{
					List[i][j] = -1;
				}
				deg[i] = 0;
				hi[i] = 0;
				spin[i] = 0;
			}

			for(i=0;i<((int)SIZE*c);i++)
			{
				for(j=0;j<2;j++)
				{
					pair[i][j] = -1;
				}
			}

			for(i=0;i<((int)SIZE*c);i++)
			{
				pair[i][0] = (int)(dsfmt_genrand_close_open(&dsfmt)*(SIZE));
				pair[i][1] = (int)(dsfmt_genrand_close_open(&dsfmt)*(SIZE));
				if(pair[i][0]==pair[i][1])
				{
					i--;
					continue;
				}
			}

			for(j=0;j<((int)SIZE*c);j++)
			{
				i = (int)((double)dsfmt_genrand_close_open(&dsfmt)*SIZE);
				if(i!=pair[j][0] && i!=pair[j][1])
				{
					List[i][deg[i]] = pair[j][0];
					deg[i] += 1;
					List[i][deg[i]] = pair[j][1];
					deg[i] += 1;
 				}else{
					j--;
					continue;
				}
			}

		//graphの構成完了
	
			for(i=0;i<SIZE;i++)
			{
				//spin[i] = 1;

				probability = (double)dsfmt_genrand_close_open(&dsfmt);
				if( probability <= (9.0/10.0) )
				{
					spin[i] = 1;
				}else if( probability > (9.0/10.0) )
				{
					spin[i] = -1;
				}

				//spin[i] = -1;
			}

		//初期状態の設定完了

		//以下モンテカルロ計算

			magnetization = 0.0;
		
			//以下のfor文でスピンを更新
			for(mcs=1;mcs<=MCS;mcs++) 
			{
				for(cite=0;cite<SIZE;cite++)//
				{
					i1=(int)((double)dsfmt_genrand_close_open(&dsfmt)*SIZE);
			
					hi[i1] = 0;
					for(j=0;j<deg[i1];j+=2) 
					{
						hi[i1] += spin[List[i1][j]]*spin[List[i1][j+1]];
					}

					r = 1.0/(1.0+exp(-2.0*hi[i1]/temperature));

					//if(hi[i1]>0){
					//	r = 1.0;
					//}else if(hi[i1]==0){
					//	r = 0.5;
					//}else{
					//	r = 0.0;
					//}
		
					if( (double)dsfmt_genrand_close_open(&dsfmt) < r )
					{
						spin[i1] = 1;
					}else{
						spin[i1] = -1;
					}

				}  //cite

				mcurrent = 0.0;
			
				for(i=0;i<SIZE;i++)
				{
					mcurrent += (double)spin[i];
				} //i

				mcurrent /= SIZE;  //1mcs あたりの磁化密度を計算

				magnetization += mcurrent;

				denominator = (double)mcs; //ここが標準誤差(errorbar)を求めるときのサンプル数(n)になる。

				Sample[mcs-1][sample] = magnetization/denominator;

			} //mcs

		} //sample

	/*配列Sample[i][j]が得られた*/

		get_average(Sample);			

	}

	for(i=0;i<MCS;i++)
	{
		free((void *)Sample[i]);
	}
	free((void *)Sample);

	for(i=0;i<SIZE;i++)
	{
		free((void *)List[i]);
	}
	free((void *)List);

	for(i=0;i<SIZE*c;i++)
	{
		free((void *)pair[i]);
	}
	free((void *)pair);

    free(spin);
	free(deg);
	free(hi);

   //getchar();
    return 0;
} 

void get_average(double **Sample)
{
	int i,j;
	double ave_magnetization;
	double ave_magnetization2;
	char filename[100];

	FILE *fp;

	sprintf(filename,"ThreeBodyMagnetization_size=16384_1sample_100000mcs_c=3.9_T=0.5_InitialCondition=0.9_%d.dat",STEP);
	fp = fopen(filename,"wb");

	for(i=0;i<MCS;i++)
	{
		ave_magnetization = 0.0;
		ave_magnetization2 = 0.0;

		if(i>=1&&i<=10)
		{
			for(j=0;j<SAMPLE;j++)
			{
				ave_magnetization += Sample[i][j];
				ave_magnetization2 += pow(Sample[i][j],2.0);
			} //j
			ave_magnetization /= (double)SAMPLE; //平均
			ave_magnetization2 /= (double)SAMPLE; //2乗平均
			fprintf(fp,"%d %f %f\n",i,fabs(ave_magnetization),sqrt( ave_magnetization2-pow(ave_magnetization,2.0) )/sqrt( (double)SAMPLE ));
		}else if(i>10&&i<=100)
		{
			if(i%10==0)
			{
				for(j=0;j<SAMPLE;j++)
				{
					ave_magnetization += Sample[i][j];
					ave_magnetization2 += pow(Sample[i][j],2.0);
				} //j
				ave_magnetization /= (double)SAMPLE;
				ave_magnetization2 /= (double)SAMPLE;
				fprintf(fp,"%d %f %f\n",i,fabs(ave_magnetization),sqrt( ave_magnetization2-pow(ave_magnetization,2.0) )/sqrt( (double)SAMPLE ));
			}
		}else if(i>100&&i<=1000)
		{
			if(i%100==0)
			{
				for(j=0;j<SAMPLE;j++)
				{
					ave_magnetization += Sample[i][j];
					ave_magnetization2 += pow(Sample[i][j],2.0);
				} //j
				ave_magnetization /= (double)SAMPLE;
				ave_magnetization2 /= (double)SAMPLE;
				fprintf(fp,"%d %f %f\n",i,fabs(ave_magnetization),sqrt( ave_magnetization2-pow(ave_magnetization,2.0) )/sqrt( (double)SAMPLE ));
			}
		}else if(i>1000&&i<=10000)
		{
			if(i%1000==0)
			{
				for(j=0;j<SAMPLE;j++)
				{
					ave_magnetization += Sample[i][j];
					ave_magnetization2 += pow(Sample[i][j],2.0);
				} //j
				ave_magnetization /= (double)SAMPLE;
				ave_magnetization2 /= (double)SAMPLE;
				fprintf(fp,"%d %f %f\n",i,fabs(ave_magnetization),sqrt( ave_magnetization2-pow(ave_magnetization,2.0) )/sqrt( (double)SAMPLE ));
			}
		}else if(i>10000&&i<=100000)
		{
			if(i%10000==0)
			{
				for(j=0;j<SAMPLE;j++)
				{
					ave_magnetization += Sample[i][j];
					ave_magnetization2 += pow(Sample[i][j],2.0);
				} //j
				ave_magnetization /= (double)SAMPLE;
				ave_magnetization2 /= (double)SAMPLE;
				fprintf(fp,"%d %f %f\n",i,fabs(ave_magnetization),sqrt( ave_magnetization2-pow(ave_magnetization,2.0) )/sqrt( (double)SAMPLE ));
			}
		}else if(i>100000&&i<=1000000)
		{
			if(i%100000==0)
			{
				for(j=0;j<SAMPLE;j++)
				{
					ave_magnetization += Sample[i][j];
					ave_magnetization2 += pow(Sample[i][j],2.0);
				} //j
				ave_magnetization /= (double)SAMPLE;
				ave_magnetization2 /= (double)SAMPLE;
				fprintf(fp,"%d %f %f\n",i,fabs(ave_magnetization),sqrt( ave_magnetization2-pow(ave_magnetization,2.0) )/sqrt( (double)SAMPLE ));
			} 
		}
	
	} //i

	fclose(fp);

	return;
}