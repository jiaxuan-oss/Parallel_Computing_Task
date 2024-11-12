//////////////////////////////////////////////////////////////////////////////////////
// MatrixGenerator_bin.c
// ----------------------------------------------------------------------------------
//
// Generates a N x M  matrix with random cell values and writes into a binary file
//
//
//////////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>

#define MATRIX_CELL_MAXVAL 1000

int main()
{
	int row, col;
	int i,j;
	char matrixName[256] = {0};
	unsigned int randNum = 0;

	printf("Welcome to the random matrix file generator!\n\n");

	printf("Specify the matrix file name (e.g. MA.bin): ");
	scanf("%s", matrixName);

	printf("Specify the matrix row and col values (e.g. 10 10): ");
	scanf("%d%d", &row, &col);

	srand((unsigned int)time(NULL));

	FILE *pFile = fopen(matrixName, "wb");
	fwrite(&row, sizeof(int), 1, pFile);
	fwrite(&col, sizeof(int), 1, pFile);

	for(i = 0; i < row; i++){
		for(j = 0; j < col; j++){
			randNum = 100 + ((unsigned int)rand() % MATRIX_CELL_MAXVAL);
			fwrite(&randNum, sizeof(unsigned int), 1, pFile);
		}
	}

	fclose(pFile);
	printf("Done!\n");
	return 0;
}