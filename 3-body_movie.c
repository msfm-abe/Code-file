#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include "dSFMT.h"

#define SIZE 10000
#define size (int)((SIZE*c)*0.2)
#define MCS 1
#define SAMPLE 1//sample数（グラフの数）


dsfmt_t dsfmt;

int main(void)
{
	char filename[64];
    int mcs;
	int sample;
    int i;
	int j;
    int i1;
	int cite;
	int seed;
    double temperature;
	double probability;
    double r;
    double c;
	int **List;
	int *deg;
	int *spin;
	int *hi;
	int **pair;

	FILE *fp;

    c = 3.9;
	temperature = 0.5;

	seed = (unsigned int)time(NULL);
	dsfmt_init_gen_rand(&dsfmt, seed);

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

			//probability = (double)dsfmt_genrand_close_open(&dsfmt);
			//if( probability <= (9.0/10.0) )
			//{
			//	spin[i] = 1;
			//}else if( probability > (9.0/10.0) )
			//{
			//	spin[i] = -1;
			//}

			spin[i] = -1;
		}

	//初期状態の設定完了

	//以下モンテカルロ計算
		
		for(mcs=1;mcs<=MCS;mcs++)
		{
			for(cite=1;cite<=SIZE;cite++)
			{
					sprintf(filename,"3-body_animation_size400_%d.dat",cite);
					fp = fopen(filename,"wb");

					for(i=1;i<=SIZE;i++){
						if(i%SIZE==0){
							fprintf(fp,"\n");
						}else{
							fprintf(fp,"%d ",spin[i]);
						}
					}
		

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

				fclose(fp);

			}  //cite

		} //mcs

	} //sample

	/*配列Sample[i][j]が得られた*/			

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
