#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>
#include <mpi.h>

//combine local matrices to the global matrix
void combineLocalMatrices(unsigned long long *localMatrix, unsigned long long *globalMatrix, int row_start, int row_end, int col_start, int col_end, int colC) {
    //only loop through the necessary row and column for the local matrix
    //+= is because each of the elem in global matrix initially is 0, so x+0 = x
    for(int row = row_start; row < row_end; row++){
        for (int col = col_start; col < col_end; col++){
            // Add localMatrix to globalMatrix element by element
            globalMatrix[row * colC + col] += localMatrix[row * colC + col];
        }
    }
}


int main(int argc, char *argv[])
{
    // Variables
    int i = 0, j = 0, k = 0;
    int my_rank, size;
    int position;
    int pack_size;
    int row_start_point, row_end_point;
    int col_start_point, col_end_point;
    unsigned long long *pMatrixC;
    int rowA, rowB, rowC, colA, colB, colC, commonPoint;
    int *pMatrixA = NULL, *pMatrixB = NULL;


    /* Clock information */
    struct timespec start, end, startComp, endComp;
    double time_taken;

    clock_gettime(CLOCK_MONOTONIC, &start);
    //initialise mpi
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    

    // 1. Read Matrix A
    rowA = colA = rowB = colB = 0;
    if (my_rank == 0) {
        printf("Matrix Multiplication using 1-Dimension Arrays - Start\n\n");

        // Reading Matrix A
        printf("Reading Matrix A - Start\n");
        FILE *pFileA = fopen("MA_500x500.bin", "rb");
        if (!pFileA) {
            perror("File A opening failed");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        fread(&rowA, sizeof(int), 1, pFileA);
        fread(&colA, sizeof(int), 1, pFileA);
        pMatrixA = (int *)malloc((rowA * colA) * sizeof(int));
        if (pMatrixA == NULL) {
            perror("Memory allocation for Matrix A failed");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        for (i = 0; i < rowA; i++) {
            //reading row by row using fread
            fread(&pMatrixA[i * colA], sizeof(int), colA, pFileA);
        }
        fclose(pFileA);
        printf("Reading Matrix A - Done\n");

        // Reading Matrix B
        printf("Reading Matrix B - Start\n");
        FILE *pFileB = fopen("MB_500x500.bin", "rb");
        if (!pFileB) {
            perror("File B opening failed");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        fread(&rowB, sizeof(int), 1, pFileB);
        fread(&colB, sizeof(int), 1, pFileB);
        pMatrixB = (int *)malloc((rowB * colB) * sizeof(int));
        if (pMatrixB == NULL) {
            perror("Memory allocation for Matrix B failed");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        for (i = 0; i < rowB; i++) {
            //reading row by row
            fread(&pMatrixB[i * colB], sizeof(int), colB, pFileB);
        }
        fclose(pFileB);
        printf("Reading Matrix B - Done\n");

        // Initialize Matrix C
        rowC = rowA;
        colC = colB;
        commonPoint = colA; // Initialize commonPoint
    }
    //broadcast rows to all the processses
    MPI_Bcast(&rowA, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&colA, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rowB, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&colB, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rowC, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&colC, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&commonPoint, 1, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast commonPoint
    //assign local matrix result
    pMatrixC = (unsigned long long *)calloc((rowC * colC), sizeof(unsigned long long));
    
    //calculate the tile partition
    //using result will divide into half even number process left odd number process right
	int threadColDiv = colC / 2;
    int threadColDivRemain = colC % 2;
    int threadRowDiv = rowC / (size / 2);
    int threadRowDivRemain = rowC % (size / 2);

    if (my_rank == 0) {
        for (int temprank = 1; temprank < size; temprank++) {
            
            position = 0;

            if (temprank % 2 == 0) {
                // Even thread
                if (temprank == (size - 2)) {
                    // Last even thread
                    row_start_point = (temprank / 2) * threadRowDiv;
                    row_end_point = row_start_point + threadRowDiv + threadRowDivRemain;
                } else {
                    // Not last even thread
                    row_start_point = (temprank / 2) * threadRowDiv;
                    row_end_point = row_start_point + threadRowDiv;
                }
                col_start_point = 0;
                col_end_point = threadColDiv;
            } else {
                // Odd thread
                if (temprank == (size - 1)) {
                    // Last odd thread
                    row_start_point = (temprank / 2) * threadRowDiv;
                    row_end_point = row_start_point + threadRowDiv + threadRowDivRemain;
                } else {
                    // Not last odd thread
                    row_start_point = (temprank / 2) * threadRowDiv;
                    row_end_point = row_start_point + threadRowDiv;
                }
                col_start_point = threadColDiv;
                col_end_point = col_start_point + threadColDiv + threadColDivRemain;
            }
            //send to all processes 
            //after getting their own row_start_point, row end  point, col start and col end
            MPI_Send(&row_start_point, 1, MPI_INT, temprank, 0, MPI_COMM_WORLD);
            MPI_Send(&row_end_point, 1, MPI_INT, temprank, 1, MPI_COMM_WORLD);
            MPI_Send(&col_start_point, 1, MPI_INT, temprank, 2, MPI_COMM_WORLD);
            MPI_Send(&col_end_point, 1, MPI_INT, temprank, 3, MPI_COMM_WORLD);
            
            //create subarray to store all the necessary rows and column 
            //needed for multiplication for that process
            int subarray_size = rowC*colC;
            int *subarrayA = (int *)malloc(subarray_size * sizeof(int));
            int *subarrayB = (int *)malloc(subarray_size * sizeof(int));
            if (!subarrayA || !subarrayB) {
                perror("Subarray allocation failed");
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }

            //loop through matrix A to copy all the necessary cells
            for (i = row_start_point; i < row_end_point; i++) {
                for (k = 0; k < commonPoint; k++){
                    subarrayA[((i - row_start_point) * colA) + k] = pMatrixA[(i * colA) + k];
                }
            }
            
            //loop through matrix B to copy all the necessary cells
            for(j = col_start_point; j < col_end_point; j++){
                for (k = 0; k < commonPoint; k++) {
                    subarrayB[((j - col_start_point) * rowA) + k] = pMatrixB[(k * colB) + j];
                }
            }
            
            //pack size according to the subarray size
            MPI_Pack_size(subarray_size + (sizeof(unsigned long long) * 2), MPI_INT, MPI_COMM_WORLD, &pack_size);
            int buffer_size = pack_size;
            //initialise buffer for subarray to send to process
            char *bufferA = (char *)malloc(buffer_size);
            char *bufferB = (char *)malloc(buffer_size);
            clock_gettime(CLOCK_MONOTONIC, &startComp); 

            //pack the time and subarrays
            MPI_Pack(&startComp.tv_sec, 1, MPI_UINT64_T, bufferA, buffer_size, &position, MPI_COMM_WORLD);
            MPI_Pack(&startComp.tv_nsec, 1, MPI_UINT64_T, bufferA, buffer_size, &position, MPI_COMM_WORLD);
            MPI_Pack(subarrayA, subarray_size, MPI_INT, bufferA, buffer_size, &position, MPI_COMM_WORLD);
            
            //send the subarrays and time together to process
            MPI_Send(bufferA, position, MPI_PACKED, temprank, 4, MPI_COMM_WORLD);
            position = 0;
            MPI_Pack(subarrayB, subarray_size, MPI_INT, bufferB, buffer_size, &position, MPI_COMM_WORLD);
            MPI_Send(bufferB, position, MPI_PACKED, temprank, 5, MPI_COMM_WORLD);

        }
    } else {
        //for other processes receive from root rank
        //receive their own row/col start and end
        MPI_Recv(&row_start_point, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&row_end_point, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&col_start_point, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&col_end_point, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //preparation to receive the packed subarrays
        int totalElem = rowC*colC + (sizeof(unsigned long long) * 2);
        MPI_Pack_size(totalElem, MPI_INT, MPI_COMM_WORLD, &pack_size);
        int buffer_size = pack_size;
        char *recv_bufferA = (char *)malloc(buffer_size);
        char *recv_bufferB = (char *)malloc(buffer_size);


        if (!recv_bufferA || !recv_bufferB) {
            perror("Buffer allocation failed");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }


        //receive the subarrays using buffer
        MPI_Recv(recv_bufferA, buffer_size, MPI_PACKED, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(recv_bufferB, buffer_size, MPI_PACKED, 0, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        clock_gettime(CLOCK_MONOTONIC, &endComp); 

        pMatrixA = (int *)malloc(totalElem * sizeof(int)); // Allocate memory for pMatrixA
        pMatrixB = (int *)malloc(totalElem * sizeof(int)); // Allocate memory for pMatrixB

        if (!pMatrixA || !pMatrixB) {
            perror("Memory allocation for pMatrixA failed");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        position = 0;
        // Unpack timestamp and subarrays
        MPI_Unpack(recv_bufferA, buffer_size, &position, &startComp.tv_sec, 1, MPI_UINT64_T, MPI_COMM_WORLD);
        MPI_Unpack(recv_bufferA, buffer_size, &position, &startComp.tv_nsec, 1, MPI_UINT64_T, MPI_COMM_WORLD);
        MPI_Unpack(recv_bufferA, buffer_size, &position, pMatrixA, (rowC*colC), MPI_INT, MPI_COMM_WORLD);
        position = 0;
        MPI_Unpack(recv_bufferB, buffer_size, &position, pMatrixB, totalElem, MPI_INT, MPI_COMM_WORLD);
        time_taken = (endComp.tv_sec - startComp.tv_sec) * 1e9; 
        time_taken = (time_taken + (endComp.tv_nsec - startComp.tv_nsec)) * 1e-9; 

        //print out the timetaken
        printf("Rank %d took %lf (s) to receive matrix from Root process\n", my_rank, time_taken);
  
    }

    if(my_rank == 0) { //rank 0 follow the normal way as the matrix arrange same as before
        row_start_point = 0;
        row_end_point = row_start_point + threadRowDiv;
        col_start_point = 0;
        col_end_point = threadColDiv;

        for (i = row_start_point; i < row_end_point; i++) {
            for (j = col_start_point; j < col_end_point; j++) {
                    pMatrixC[i * colC + j] = 0; // Initialize the result element
                for (k = 0; k < commonPoint; k++) {
                    pMatrixC[(i*colC)+j] += (pMatrixA[(i*colA)+k] * pMatrixB[(k*colB)+j]);
                }
            }
        }
    } else{ //the others the other way as it is transposed

        for (i = row_start_point; i < row_end_point; i++) {
            for (j = col_start_point; j < col_end_point; j++) {
                pMatrixC[i * colC + j] = 0; // Initialize the result element
                for (k = 0; k < commonPoint; k++) {
                    // Matrix A row remains the same
                    int elemA = pMatrixA[((i - row_start_point) * commonPoint) + k];

                    // Matrix B column is treated as a row
                    int elemB = pMatrixB[((j - col_start_point) * commonPoint) + k];

                    // Perform the multiplication
                    pMatrixC[(i * colC) + j] += elemA * elemB;
                }
            }
        }
    }
    
    unsigned long long *globalMatrixC = NULL;
    if(my_rank != 0){
        //send local pMatrixC back to root
        MPI_Pack_size(rowC*colC, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD, &pack_size);
        int buffer_size = pack_size;
        char *bufferC = (char *)malloc(buffer_size);
        position = 0;
        MPI_Pack(pMatrixC, rowC * colC, MPI_UNSIGNED_LONG_LONG, bufferC, buffer_size, &position, MPI_COMM_WORLD);
        MPI_Send(bufferC, position, MPI_PACKED, 0, my_rank, MPI_COMM_WORLD);
        MPI_Send(&row_start_point, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&row_end_point, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Send(&col_start_point, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
        MPI_Send(&col_end_point, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
            
    } 

    if (my_rank == 0){
        globalMatrixC = (unsigned long long *)calloc((rowC * colC), sizeof(unsigned long long));
        if (!globalMatrixC) {
            perror("Memory allocation for globalMatrixC failed");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        //set all values to 0
        for(i = 0; i < rowC* colC; i++){
            globalMatrixC[i] = 0;
        }
        //combine the rank 0 matrixC to globalmatrixC
        combineLocalMatrices(pMatrixC, globalMatrixC, row_start_point, row_end_point, col_start_point, col_end_point, colC);
       
        //loop through each rank and receive their local matrixC from each rank
        //then combine it to globalmatrixC
        for (i = 1; i < size; i++){


            MPI_Pack_size(rowC*colC, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD, &pack_size);
            int buffer_size = pack_size;
            char *recv_bufferC = (char *)malloc(buffer_size);

            //receive bufferC from each rank
            MPI_Recv(recv_bufferC, buffer_size, MPI_PACKED, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            position = 0;
            MPI_Unpack(recv_bufferC, buffer_size, &position, pMatrixC, rowC*colC, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
            MPI_Recv(&row_start_point, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&row_end_point, 1, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&col_start_point, 1, MPI_INT, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&col_end_point, 1, MPI_INT, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            
            combineLocalMatrices(pMatrixC, globalMatrixC, row_start_point, row_end_point, col_start_point, col_end_point, colC);
            
        }
    }
    //checking each rank's row/col start and end
    // printf("rank: %d, row_start: %d, row_end: %d, col_start: %d, col_end: %d", my_rank, row_start_point, row_end_point, col_start_point, col_end_point);
    // fflush(stdout);
    
    //waiting for each rank to finish
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    if (my_rank == 0) {
        // Write result to file
        printf("Writing Matrix C - Start\n");
        FILE *pFileC = fopen("MC_500x500.bin", "wb");
        fwrite(&rowC, sizeof(int), 1, pFileC);
        fwrite(&colC, sizeof(int), 1, pFileC);
        //write row by row
        for (i = 0; i < rowC; i++) {
            fwrite(&globalMatrixC[i * colC], sizeof(unsigned long long), colC, pFileC);
        }

    
        fclose(pFileC);
        printf("Writing Matrix C - Done\n");

        // Time calculation of overall program
        clock_gettime(CLOCK_MONOTONIC, &end);
        time_taken = (end.tv_sec - start.tv_sec) * 1e9;
        time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
        printf("Elapsed Time: %f seconds\n", time_taken);
    }

    return 0;
}
