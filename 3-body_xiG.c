#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"dSFMT.h"
 
#define SIZE 100
#define EDGE SIZE*c
#define size (int)((SIZE*c)*0.2)
#define MCS 5000
#define SAMPLE 2
#define GRAPH 2
#define InitialMag 1.0
#define TEMPERATURE 0.5
 
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
    int mcs;
    int sample;
    int graph;
    double probability;
    int hi;
    double mcurrent;
    double **memory;
    double **tr_average;
    double ave_magnetization;
    double ave_magnetization2;
    double xiG;
 
    FILE *fp;
    fp = fopen("xiG_test.dat","wb");
     
    seed = (unsigned int)time(NULL);
    dsfmt_init_gen_rand(&dsfmt,seed);
 
    c = 1.0;
 
 
    spin = (int *)malloc(sizeof(int)*SIZE);
 
    memory = (double **)malloc(sizeof(double *)*SAMPLE);
    for(i=0;i<SAMPLE;i++){
        memory[i] = (double *)malloc(sizeof(double)*MCS);
    }
 
    tr_average = (double **)malloc(sizeof(double *)*GRAPH);
    for(i=0;i<GRAPH;i++){
        tr_average[i] = (double *)malloc(sizeof(double)*MCS);
    }
 
    for(i=0;i<GRAPH;i++){
        for(j=0;j<MCS;j++){
            tr_average[i][j] = 0.0;
        }
    }
 
 
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
 
        //for(i=0;i<EDGE;i++){
        //  printf("pair[%d][0]=%d pair[%d][1]=%d\n",i,pair[i][0],i,pair[i][1]);
        //}
        //printf("\n");
        //for(i=0;i<SIZE;i++){
        //  printf("List[%d] %d ",i,deg[i]);
        //  for(j=0;j<size;j++){
        //      printf("%d ",List[i][j]);
        //  }
        //  printf("\n");
        //}
        //      
 
        for(i=0;i<SAMPLE;i++){
            for(j=0;j<MCS;j++){
                memory[i][j] = 0.0;
            }
        }
         
        //for(i=0;i<SAMPLE;i++){
        //  for(j=0;j<MCS;j++){
        //      printf("memory[%d][%d]=%f\n",i,j,memory[i][j]);
        //  }
        //}
 
        for(sample=0;sample<SAMPLE;sample++){
 
            for(i=0;i<SIZE;i++){
                if( (double)dsfmt_genrand_close_open(&dsfmt)<=(InitialMag+1.0)/2.0 ){
                    spin[i] = 1;
                }else{
                    spin[i] = -1;
                }
            }
 
            //for(i=0;i<SIZE;i++){
            //  printf("spin[%d]=%d\n",i,spin[i]);
            //}
            //printf("\n");
 
            magnetization = 0.0;
 
            for(mcs=1;mcs<=MCS;mcs++){
 
                for(step=0;step<SIZE;step++){
                    choice = (int)(dsfmt_genrand_close_open(&dsfmt)*(SIZE));
             
                    hi = 0;
                    for(i=0;i<deg[choice];i+=2){
                        hi += spin[List[choice][i]]*spin[List[choice][i+1]];
                    }
 
                    probability = 1.0/( 1.0+exp(-2.0*hi/TEMPERATURE) );
 
                    if( (double)dsfmt_genrand_close_open(&dsfmt)<=probability ){
                        spin[choice] = 1;
                    }else{
                        spin[choice] = -1;
                    }
                }//step
 
                //for(i=0;i<SIZE;i++){
                //  printf("spin[%d]=%d\n",i,spin[i]);
                //}
 
                mcurrent = 0.0;
                for(i=0;i<SIZE;i++){
                    mcurrent += (double)spin[i];
                }
                //printf("mcurrent=%f\n",mcurrent);
                mcurrent /= (double)SIZE;
                //printf("mcurrent/SIZE=%f\n",mcurrent);
 
                magnetization += mcurrent;
 
                memory[sample][mcs-1] = magnetization/mcs;
                //printf("memory[%d][%d]=%f\n magnetization=%f\n",sample,mcs-1,memory[sample][mcs-1],magnetization);
                //printf("\n");
             
            }//mcs
            //printf("\n");
        }//sample
 
        for(i=0;i<MCS;i++){
            for(j=0;j<SAMPLE;j++){
                tr_average[graph][i] += memory[j][i];
            }
            tr_average[graph][i] /= SAMPLE;
            //printf("tr_average[%d][%d]=%f ",graph,i,tr_average[graph][i]);
        }
        //printf("\n");
     
    }//graph
 
    ////derive a xiG
    for(i=0;i<MCS;i++){
 
        //if( (i%10000)==9999 ){
            ave_magnetization = 0.0;
            ave_magnetization2 = 0.0;
            xiG = 0.0;
 
            for(graph=0;graph<GRAPH;graph++){
                ave_magnetization += tr_average[graph][i];
                ave_magnetization2 += pow(tr_average[graph][i],2);
                //printf("ave_magnetization=%f\n",ave_magnetization);
                //printf("ave_magnetization2=%f\n",ave_magnetization2);
            }
            ave_magnetization /= GRAPH;
            ave_magnetization2 /= GRAPH;
 
            xiG = SIZE*(ave_magnetization2-pow(ave_magnetization,2));
            //printf("xiG=%f\n",xiG);
 
            fprintf(fp,"%d %f %f\n",i+1,fabs(xiG),sqrt( (ave_magnetization2-pow(ave_magnetization,2)) )/sqrt( (double)GRAPH) );
 
        //}
    }
 
    fclose(fp);
 
    free((void *)spin);
    for(i=0;i<SIZE;i++){
        free((void *)List[i]);
    }
    free((void *)List);
    free((void *)deg);
    for(i=0;i<EDGE;i++){
        free((void *)pair[i]);
    }
    free((void *)pair);
    for(i=0;i<SAMPLE;i++){
        free((void *)memory[i]);
    }
    free((void *)memory);
    for(i=0;i<GRAPH;i++){
        free((void *)tr_average[i]);
    }
    free((void *)tr_average);
 
 
    //getchar();
    return 0;
 
}