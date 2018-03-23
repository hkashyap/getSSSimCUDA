/*
* getSSSim version used in the Scientific reports paper
* Author: Hirak Kashyap
* Date: 01/24/2016
* Update Makefile to link to this source file
 */
#include <sys/time.h>
#include <stdio.h>
#include <cuda_runtime.h>
#include <string.h>
#include <math.h>

#define nrows 9084
#define ncols 28

#define sBuffer 1024

float h_data[nrows][ncols];
float h_score[nrows][nrows];


/**
 * CUDA Kernel Device code
 *
 */

__global__ void getSSSim(const float *data, float *score)
{

	int tid_x = blockIdx.x*blockDim.x + threadIdx.x;
	int tid_y = blockIdx.y*blockDim.y + threadIdx.y;
	
	int s_grid_width = gridDim.x * blockDim.x;
	int s_index = s_grid_width * tid_y + tid_x;

	int d_grid_width = gridDim.z * blockDim.z;
	int dx_index = d_grid_width * tid_x;
	int dy_index = d_grid_width * tid_y;
			
	float *gx_d = (float*)((char*) data + dx_index * sizeof(float));
	float *gy_d = (float*)((char*) data + dy_index * sizeof(float));

	//__shared__ float sum;
	//float sum = 0;
	float score_s = 0;
	float score_l[27];
	float gx[27];
	float gy[27];
	float gx_2min1 = gx_d[1] - gx_d[0];
	float gy_2min1 = gy_d[1] - gy_d[0];

	//0-26
	for(int i=0; i<ncols-1; i++)
	{
		gx[i] = (gx_d[i+1] - gx_d[i])/gx_2min1;
		gy[i] = (gy_d[i+1] - gy_d[i])/gy_2min1;
		
	}
	__syncthreads();

	//0 and 26
	score_l[0]=(gx[0]+gx[1]+gy[0]+gy[1])/4;
    	score_l[ncols-2]=(gx[ncols-3]+gx[ncols-2]+gy[ncols-3]+gy[ncols-2])/4;

    	//1-25
    	for(int i=1; i<ncols-2; i++)
    	{
        	score_l[i]=(gx[i-1]+gx[i]+gx[i+1]+gy[i-1]+gy[i]+gy[i+1])/6;
    	}
//	__syncthreads();
    	int n_diff=0;
    	for(int i=0; i<ncols-1; i++)
    	{
        	score_l[i] = 2 * fmaxf(fabsf(score_l[i]-gx[i]),fabsf(score_l[i]-gy[i]));
        	if (score_l[i]>=0.00001)
		{
        	      n_diff++;
        	      score_s = score_s + fabsf(gx[i]-gy[i])/score_l[i];
		}
    	}
	__syncthreads();
    //	score[s_index] = 1 - (score_s/ncols);
    	score[s_index] = 1 - (score_s/n_diff);
	__syncthreads();
}

void populateArrays(){
    char buf[sBuffer];
    FILE *fp;
    fp = fopen("ADNormal_onlyExp.csv", "r");
    for(int i=0;i<nrows;i++)
    {
	fgets(buf, sizeof(buf), fp);
	char *tok = strtok(buf,",");
	h_data[i][0] = atof(tok);
	for(int j=1;j<ncols;j++)
	{
		tok = strtok(NULL, ",");
		h_data[i][j] = atof(tok);
	}
    }
}

void print_cal()
{
	float maxdiff = 0;
	//float score[nrows][nrows];
	float score;
	for(int i=0; i< nrows; i++)
	{
		
		//printf("\n");
		for(int j=0; j<nrows; j++)
		{
			float score_s = 0;
			float score_l[27];
			float gx_2min1 = h_data[j][1] - h_data[j][0];
			float gy_2min1 = h_data[i][1] - h_data[i][0];
			//float gy_2min1 = gy[1] - gy[0];
			float gx[27], gy[27];

			for(int k=0; k<ncols-1; k++)
			{
				gx[k] = (h_data[j][k+1] - h_data[j][k])/gx_2min1;
				gy[k] = (h_data[i][k+1] - h_data[i][k])/gy_2min1;
			}
			//__syncthreads();

			//0 and 26
			score_l[0]=(gx[0]+gx[1]+gy[0]+gy[1])/4;
    			score_l[ncols-2]=(gx[ncols-3]+gx[ncols-2]+gy[ncols-3]+gy[ncols-2])/4;

    			//1-25
    			for(int k=1; k<ncols-2; k++)
    			{
        			score_l[k]=(gx[k-1]+gx[k]+gx[k+1]+gy[k-1]+gy[k]+gy[k+1])/6;
    			}
			//__syncthreads();
    			int n_diff=0;
    			for(int k=0; k<ncols-1; k++)
    			{
        			score_l[k] = 2 * fmaxf(fabsf(score_l[k]-gx[k]),fabsf(score_l[k]-gy[k]));
        			if (score_l[k]>=0.00001)
				{
        			      n_diff++;
        			      score_s = score_s + fabsf(gx[k]-gy[k])/score_l[k];
				}
				//__syncthreads();
    			}
    			//score = 1 - (score_s/(ncols-1));
    			score = 1 - (score_s/n_diff);
			
			//printf("%f - %f\t",score, h_score[i][j]);
			if (abs(score - h_score[i][j])>maxdiff)
				maxdiff = abs(score - h_score[i][j]); 
	/*	*/	
		}
	}
	printf("\nMax diff %f", maxdiff);
}

void saveResults()
{
	FILE *fp;
	
	fp =fopen("results_nonan.csv","w");
	for(int i=0; i<nrows; i++)
	{
		float score;
		for(int j=0; j<nrows-1; j++)
		{
			score = h_score[i][j];
			if(isnan(score)){
				score = 0.0;
			}
			fprintf(fp, "%f,", score);
		}	
		score = h_score[i][nrows-1];
		if(isnan(score)){
			score = 0.0;
		}
		fprintf(fp, "%f\n", score);
	}
	fclose(fp);
}

/**
 * Host main routine
 */
int
main(void)
{
    printf("[ getSSSim of all genes ]");

    // Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;
    timeval t0, t1, t0_cpu, t1_cpu;

    //size_t height = nrows;
    //size_t w_size = ncols * sizeof(float);
    //size_t h_size = nrows * sizeof(float);

    // Allocate the host output score matrix

    // Pitches for data and score
    //size_t pitch_data, pitch_score;

    populateArrays();

    //starts actual execution
    gettimeofday(&t0, 0);

    // Allocate the device input data matrix
    float *d_data = NULL;
    //err = cudaMallocPitch((void **)&d_data, &pitch_data, w_size, nrows);
    err = cudaMalloc((void **)&d_data, nrows*ncols*sizeof(float));

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vector d_data (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Allocate the device output score matrix
    float *d_score = NULL;
    err = cudaMalloc((void **)&d_score, nrows*nrows*sizeof(float));
    //err = cudaMallocPitch((void **)&d_score, &pitch_score, h_size, nrows);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vector d_score (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    //Copy input data from the host memory to the CUDA device
    //err = cudaMemcpy(d_g1, h_g1, size, cudaMemcpyHostToDevice);
    err = cudaMemcpy(d_data, h_data, nrows*ncols*sizeof(float), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy vector g1 from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    //launch the getSSSim kernel 
    dim3 blockSize;
    blockSize.x = 4;
    blockSize.y = 4;
    blockSize.z = ncols;

    dim3 gridSize;
    gridSize.x = nrows/blockSize.x;
    gridSize.y = nrows/blockSize.y;



    getSSSim<<<gridSize, blockSize>>>(d_data, d_score);
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to launch vectorAdd kernel (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy the device result vector in device memory to the host result vector
    // in host memory.

    //err = cudaMemcpy2D(h_score, h_size, d_score, pitch_score, h_size, nrows, cudaMemcpyDeviceToHost);
    err = cudaMemcpy(h_score, d_score, nrows*nrows*sizeof(float), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy vector g1 from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }


    gettimeofday(&t1, 0);
    long long elapsed = (t1.tv_sec-t0.tv_sec)*1000000LL + t1.tv_usec-t0.tv_usec;
    printf("\nTime elapsed GPU (microsecond):%lld\n", elapsed);
    
    printf("\nCPU\n");
    gettimeofday(&t0_cpu, 0);
    print_cal();
    gettimeofday(&t1_cpu, 0);

    
    long long elapsed_cpu = (t1_cpu.tv_sec-t0_cpu.tv_sec)*1000000LL + t1_cpu.tv_usec-t0_cpu.tv_usec;
    printf("\nTime elapsed CPU (microsecond):%lld\n", elapsed_cpu);

    //Save results into a csv file
    saveResults();

    // Free device global memory
    err = cudaFree(d_data);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to free device vector data (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    err = cudaFree(d_score);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to free device vector score (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Free host memory
    //free(h_data);
    //free(h_score);
    //free(h_score);

    // Reset the device and exit
    // cudaDeviceReset causes the driver to clean up all state. While
    // not mandatory in normal operation, it is good practice.  It is also
    // needed to ensure correct operation when the application is being
    // profiled. Calling cudaDeviceReset causes all profile data to be
    // flushed before the application exits
    err = cudaDeviceReset();

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to deinitialize the device! error=%s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    printf("Done\n");
    return 0;
}

