/*
* getSSSim CPU-only mode 
* Author: Hirak J Kashyap
* Date: March 08, 2018
*/
#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define nrows 907
#define ncols 55

#define sBuffer 10240

float h_data[nrows][ncols];
//float h_score[nrows][nrows];
float c_score[nrows][nrows];
float c_row_sum[nrows];

void populateArrays(){
    char buf[sBuffer];
    FILE *fp;

  
    fp = fopen("GSE72962_PD_miRNA_SSSim.csv", "r");
    //fp = fopen("ADNormal_pp_1.csv", "r");
    for(int i=0;i<nrows;i++)
    {
		printf("\n");
		fgets(buf, sizeof(buf), fp);
		char *tok = strtok(buf,",");
		h_data[i][0] = atof(tok);
		printf("%f\t",h_data[i][0]);
		for(int j=1;j<ncols;j++)
		{
			tok = strtok(NULL, ",");
			h_data[i][j] = atof(tok);
			if (j<10){
				printf("%f\t",h_data[i][j]);
			}
		}
		printf("\n");
    }
}

void print_cal()
{
	for(int i=0; i< nrows; i++)
	{
		float row_sum = 0;
		for(int j=0; j<nrows; j++)
		{
			float score_s = 0;
			float score_l[ncols-1];
			float gx_2min1 = h_data[j][1] - h_data[j][0];
			float gy_2min1 = h_data[i][1] - h_data[i][0];
			float gx[ncols-1], gy[ncols-1];

			for(int k=0; k<ncols-1; k++)
			{
				gx[k] = (h_data[j][k+1] - h_data[j][k])/gx_2min1;
				gy[k] = (h_data[i][k+1] - h_data[i][k])/gy_2min1;
			}
			score_l[0]=(gx[0]+gx[1]+gy[0]+gy[1])/4;
    			score_l[ncols-2]=(gx[ncols-3]+gx[ncols-2]+gy[ncols-3]+gy[ncols-2])/4;

    			for(int k=1; k<ncols-2; k++)
    			{
        			score_l[k]=(gx[k-1]+gx[k]+gx[k+1]+gy[k-1]+gy[k]+gy[k+1])/6;
    			}
    			float n_diff=0.0;
    			for(int k=0; k<ncols-1; k++)
    			{
        			score_l[k] = 2.0 * fmaxf(fabsf(score_l[k]-gx[k]),fabsf(score_l[k]-gy[k]));
        			if (score_l[k]>=0.0000001)
				{
        			      n_diff = n_diff + 1.0;
        			      score_s = score_s + fabsf(gx[k]-gy[k])/score_l[k];
				}
    			}
    			c_score[i][j] = 1.0 - (score_s/n_diff);
			
			if(i==j){
				c_score[i][j] = 1.0;
			}	
			
			if(isnan(c_score[i][j])){
				c_score[i][j] = 0.0;
			}
			
			row_sum = row_sum + c_score[i][j];
		}
		c_row_sum[i] = row_sum;
	}
	
}

void saveResults()
{
	FILE *fp;
	
	fp =fopen("results_GSE72962_PD_miRNA_SSSim.csv","w");
	for(int i=0; i<nrows; i++)
	{
		float score;
		for(int j=0; j<nrows-1; j++)
		{
			score = c_score[i][j];
			fprintf(fp, "%f,", score);
		}	
		score = c_score[i][nrows-1];
		if(isnan(score)){
			score = 0.0;
		}
		fprintf(fp, "%f\n", score);
	}
	fclose(fp);

	fp =fopen("results_rowsum_GSE72962_PD_miRNA_SSSim.csv","w");
	for(int i=0; i<nrows; i++)
	{
		fprintf(fp, "%f\n", c_row_sum[i]);
	}
	fclose(fp);
}

/**
 * Host main routine
 */
int main(void)
{
    printf("[ getSSSim of all genes ]");

    // Error code to check return values for CUDA calls
    //cudaError_t err = cudaSuccess;
    timeval t0_cpu, t1_cpu;


    populateArrays();

    
    printf("\nCPU\n");
    gettimeofday(&t0_cpu, 0);
    print_cal();
    gettimeofday(&t1_cpu, 0);

    
    long long elapsed_cpu = (t1_cpu.tv_sec-t0_cpu.tv_sec)*1000000LL + t1_cpu.tv_usec-t0_cpu.tv_usec;
    printf("\nTime elapsed CPU (microsecond):%lld\n", elapsed_cpu);

    //Save results into a csv file
    saveResults();

    printf("Done\n");
    return 0;
}
