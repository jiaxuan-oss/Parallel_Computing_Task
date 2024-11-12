//////////////////////////////////////////////////////////////////////////////////////
// MatrixMul_1D_bin.c
// ----------------------------------------------------------------------------------
//
// Multiplies two matrices and writes the resultant multiplication into a binary file.
//
//////////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>
#include <pthread.h>
#define NUM_THREADS 8

int *pMatrixA, *pMatrixB;
unsigned long long *pMatrixC;
int rowA, rowB, rowC, colA, colB, colC, commonPoint;

void *ThreadFunc2(void *pArg){
    // Assume the matrix arrays, row, column and commonPoint variables are global variables and initialised in main function. 
    int i, j, k;
    int my_rank = *((int*)pArg);
    
	//calculate the tile partition
    //using result will divide into half even number process left odd number process right
    int row_start_point, row_end_point;
    int col_start_point, col_end_point;
    int threadColDiv = 1;
    int threadColDivRemain = 0;
    int threadRowDiv = 1;
    int threadRowDivRemain = 0;
    
    threadColDiv = colC / 2;
    threadColDivRemain = colC % 2;
    threadRowDiv = rowC / (NUM_THREADS / 2);
    threadRowDivRemain = rowC % (NUM_THREADS / 2);

    if(my_rank % 2 == 0){
        // Even thread
        if(my_rank == (NUM_THREADS - 2)){
            // Last even thread
            row_start_point = (my_rank / 2) * threadRowDiv;
            row_end_point = row_start_point + threadRowDiv + threadRowDivRemain;
        }else{
            // Not last even thread
            row_start_point = (my_rank / 2) * threadRowDiv;
            row_end_point = row_start_point + threadRowDiv;
        }
        col_start_point = 0;
        col_end_point = threadColDiv;
    }else{
        // Odd thread
        if(my_rank == (NUM_THREADS - 1)){
            // Last odd thread
            row_start_point = (my_rank / 2) * threadRowDiv;
            row_end_point = row_start_point + threadRowDiv + threadRowDivRemain;
        }else{
            // Not last odd thread
            row_start_point = (my_rank / 2) * threadRowDiv;
            row_end_point = row_start_point + threadRowDiv;
        }
        col_start_point = threadColDiv;
        col_end_point = col_start_point + threadColDiv + threadColDivRemain;
    }
    
    // Matrix multiplication
    for(i = row_start_point; i < row_end_point; i++){
        for(j = col_start_point; j < col_end_point; j++){
            for(k = 0; k < commonPoint; k++){
                pMatrixC[(i*colC)+j] += (pMatrixA[(i*colA)+k] * pMatrixB[(k*colB)+j]);
            }
        }
    }
    return NULL;
}

int main()
{
	// Variables
	int i = 0, j = 0, k = 0;

	/* Clock information */
	struct timespec start, end; 
	double time_taken;
	pthread_mutex_t mutex;
	pthread_mutex_init(&mutex, NULL);
    pthread_t threads[NUM_THREADS];
    int thread_ids[NUM_THREADS];
	clock_gettime(CLOCK_MONOTONIC, &start); 
	
	// 1. Read Matrix A
	rowA = 0;
	colA = 0;

	printf("Matrix Multiplication using 1-Dimension Arrays - Start\n\n");

	printf("Reading Matrix A - Start\n");

	FILE *pFileA = fopen("MA_500x500.bin", "rb");
	fread(&rowA, sizeof(int), 1, pFileA); 
	fread(&colA, sizeof(int), 1, pFileA); 

	pMatrixA = (int*)malloc((rowA*colA) * sizeof(int));
	for(i = 0; i < rowA; i++){
		fread(&pMatrixA[i*colA], sizeof(int), colA, pFileA);
	}
	fclose(pFileA);

	printf("Reading Matrix A - Done\n");

	// 2. Read Matrix B
	rowB = 0;
	colB = 0;

	printf("Reading Matrix B - Start\n");
	
	FILE *pFileB = fopen("MB_500x500.bin", "rb");
	fread(&rowB, sizeof(int), 1, pFileB); 
	fread(&colB, sizeof(int), 1, pFileB); 

	pMatrixB = (int*)malloc((rowB*colB) * sizeof(int));
	for(i = 0; i < rowB; i++){
		fread(&pMatrixB[i*colB], sizeof(int), colB, pFileB); 
	}
	fclose(pFileB);

	printf("Reading Matrix B - Done\n");

	// 3. Perform matrix multiplication 
	printf("Matrix Multiplication - Start\n");

	rowC = rowA;
	colC = colB;
	pMatrixC = (unsigned long long*)calloc((rowC*colC), sizeof(unsigned long long));
    commonPoint = colC;

	int rows_per_thread = rowC/NUM_THREADS;
	int cols_per_thread = colC/NUM_THREADS;

	for (i = 0; i < NUM_THREADS; i++){
		//creating threads to compute
        thread_ids[i] = i;
		pthread_create(&threads[i], NULL, ThreadFunc2, &thread_ids[i]);
	}

	for (i = 0; i< NUM_THREADS; i++){
		//join them tgt after computing
		pthread_join(threads[i], NULL);
	}

	

	printf("Matrix Multiplication - Done\n");

	// 4. Write resuls to a new file
	printf("Write Resultant Matrix C to File - Start\n");

	FILE *pFileC = fopen("MC_500x500.bin", "wb"); 
	fwrite(&rowC, sizeof(int), 1, pFileC); 
	fwrite(&colC, sizeof(int), 1, pFileC); 
	for(i = 0; i < rowC; i++){
		fwrite(&pMatrixC[i*colC], sizeof(unsigned long long), colC, pFileC); // use fwrite to write one row of columns to the file
	}



	fclose(pFileC);
	printf("Write Resultant Matrix C to File - Done\n");
	//time taken for them to complete overall program
	clock_gettime(CLOCK_MONOTONIC, &end); 
	time_taken = (end.tv_sec - start.tv_sec) * 1e9; 
	time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9; 
	printf("Overall time (Including read, multiplication and write)(s): %lf\n", time_taken);	// ts

	// Clean up
	free(pMatrixA);
	free(pMatrixB);
	free(pMatrixC);

	printf("Matrix Multiplication using 1-Dimension Arrays - Done\n");
	return 0;
}