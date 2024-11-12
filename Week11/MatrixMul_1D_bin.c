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



int main()
{
	// Variables
	int i = 0, j = 0, k = 0;

	/* Clock information */
	struct timespec start, end, startComp, endComp; 
	double time_taken;

	clock_gettime(CLOCK_MONOTONIC, &start); 
	
	// 1. Read Matrix A
	int rowA = 0, colA = 0;

	printf("Matrix Multiplication using 1-Dimension Arrays - Start\n\n");

	printf("Reading Matrix A - Start\n");

	FILE *pFileA = fopen("MA_500x500.bin", "rb");
	fread(&rowA, sizeof(int), 1, pFileA); 
	fread(&colA, sizeof(int), 1, pFileA); 

	int *pMatrixA = (int*)malloc((rowA*colA) * sizeof(int));
	for(i = 0; i < rowA; i++){
		fread(&pMatrixA[i*colA], sizeof(int), colA, pFileA);
	}
	fclose(pFileA);

	printf("Reading Matrix A - Done\n");

	// 2. Read Matrix B
	int rowB = 0, colB = 0;

	printf("Reading Matrix B - Start\n");
	
	FILE *pFileB = fopen("MB_500x500.bin", "rb");
	fread(&rowB, sizeof(int), 1, pFileB); 
	fread(&colB, sizeof(int), 1, pFileB); 

	int *pMatrixB = (int*)malloc((rowB*colB) * sizeof(int));
	for(i = 0; i < rowB; i++){
		fread(&pMatrixB[i*colB], sizeof(int), colB, pFileB); 
	}
	fclose(pFileB);

	printf("Reading Matrix B - Done\n");

	// 3. Perform matrix multiplication 
	printf("Matrix Multiplication - Start\n");

	int rowC = rowA, colC = colB;
	unsigned long long *pMatrixC = (unsigned long long*)calloc((rowC*colC), sizeof(unsigned long long));

	int commonPoint = colA;

	clock_gettime(CLOCK_MONOTONIC, &startComp); 

	for(i = 0; i < rowC; i++){
		for(j = 0; j < colC; j++){
			for(k = 0; k < commonPoint; k++){
				pMatrixC[(i*colC)+j] += (pMatrixA[(i*colA)+k] * pMatrixB[(k*colB)+j]);
			}
		}
	}

	clock_gettime(CLOCK_MONOTONIC, &endComp);

	printf("Matrix Multiplication - Done\n");

	// 4. Write resuls to a new file
	printf("Write Resultant Matrix C to File - Start\n");

	FILE *pFileC = fopen("MC_500x500.bin", "wb"); 
	fwrite(&rowC, sizeof(int), 1, pFileC); 
	fwrite(&colC, sizeof(int), 1, pFileC); 
	for(i = 0; i < rowC; i++){
		fwrite(&pMatrixC[i*colC], sizeof(unsigned long long), colC, pFileC); // use fwrite to write one row of columns to the file
		printf("%lld ",pMatrixC[i * colC] );

	}

	

	fclose(pFileC);
	printf("Write Resultant Matrix C to File - Done\n");

	clock_gettime(CLOCK_MONOTONIC, &end); 
	time_taken = (end.tv_sec - start.tv_sec) * 1e9; 
	time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9; 
	printf("Overall time (Including read, multiplication and write)(s): %lf\n", time_taken);	// ts
	
	time_taken = (endComp.tv_sec - startComp.tv_sec) * 1e9; 
	time_taken = (time_taken + (endComp.tv_nsec - startComp.tv_nsec)) * 1e-9; 
	printf("Overall computation time(s): %lf\n", time_taken);	// ts
	

	// Clean up
	free(pMatrixA);
	free(pMatrixB);
	free(pMatrixC);

	printf("Matrix Multiplication using 1-Dimension Arrays - Done\n");
	return 0;
}