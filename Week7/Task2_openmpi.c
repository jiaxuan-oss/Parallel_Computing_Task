#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <mpi.h>

//function prototype
bool isprime(int number);
int compare(const void *a, const void *b);

int main (int argc, char *argv[]) 
{   
    //function declaration and initialisation
    int *arr = NULL;
    struct timespec start, end, startComp, endComp;
    double time_taken;
    int n = 10000000, rank, size, localcounter = 0;
    int totalCount = 0;
    int counter = 0;
    int *counts, *displacement;
    int *localarr = NULL;

    //Initialise the MPI environment and processes
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get current rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); //get all the processes num

    //only root process
    if (rank == 0){
        printf("Compute:\n");

        //counts for gathering all the count of local prime result in an array
        counts = (int *)malloc(size * sizeof(int));

        //displacement for gathering the displacement of the each result array from each processes
        displacement = (int *)malloc(size * sizeof(int));   
             
    }
    else{
        counts = NULL;
        displacement = NULL;
    }

    //localarr to store local prime number 
    localarr = (int *)malloc(n * sizeof(int));
    
    // Get current clock time to time the computational time
    clock_gettime(CLOCK_MONOTONIC, &start);
    
    //using round robin work distribution
    //use mod to distribute task -> more evenly distributed
    for (int i = 2; i < n; i++){
        if (i % size == rank){
            if (isprime(i)){
                //store in each local array of processes
                localarr[localcounter] = i;
                localcounter++;
            }
        }
    }
    
    //wait for other processes to done before proceeding
    MPI_Barrier(MPI_COMM_WORLD);

    //Gather all the local count to an array eg {23,43,15,4}
    MPI_Gather(&localcounter, 1,MPI_INT, counts, 1,MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0){
        displacement[0] = 0;
        totalCount = counts[0];
        for (int i = 1 ;i < size; i++){
            //add all the counts from each processes to total count
            totalCount += counts[i];
            //calculate the displacement of each processes 
            //using the counts to know each array have how many values
            //and then set the displacement
            //for the program to know where to start
            displacement[i] = displacement[i - 1] + counts[i - 1];
        }
    }
    if (rank == 0){
        arr = (int *)malloc(totalCount * sizeof(int));
    }
    //gather all local array into a single array using counts and displacement
    MPI_Gatherv(localarr, localcounter, MPI_INT, arr, counts, displacement, MPI_INT, 0, MPI_COMM_WORLD);
    //sort the prime array using compare function and qsort 
    qsort(arr, totalCount, sizeof(int), compare);

    //computational time end
    clock_gettime(CLOCK_MONOTONIC, &end);
    if (rank == 0){
        //print all the output
        for(int i = 0; i < totalCount; i++){
            printf("%d ", arr[i]);
        }
    }


    //get the time 
    if(rank == 0){
    time_taken = (end.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("Overall time: %f sec \n", time_taken);
    }
    MPI_Finalize();
    return 0;
}

// Comparison function for qsort
int compare(const void *a, const void *b) {
    //if it is positive number means a is larger
    //if is negative number means b is larger
    //deferencing to access the integer value
    return (*(int *)a - *(int *)b);
}



bool isprime(int number){

    if (number <= 1){
        return false;
    }

    if (number == 2 || number == 3){
        return true;
    }

    if (number % 2 == 0){
        return false;
    }
    
    if (number % 3 == 0){
        return false;
    }

    //loop until sqrt of n as p*q = n so loop until one of them is sufficient
    for (int i = 5; i <= sqrt(number); i++ ){
        if (number % i == 0){
            return false;
        }
    }
    return true;
}