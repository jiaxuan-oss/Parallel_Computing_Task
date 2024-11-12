#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char* argv[])
{
    //variable declaration and initialisation
    long N = -1;
    int i, my_rank, p;
    double sum = 0.0, total_sum = 0.0; //local sum for processor and total sum when all processor sum tgt
    double piVal; //pi value
    struct timespec start, end;
    double time_taken;
    int startIndex, endIndex, localN; //variale for processor

    //initialise MPI and get current processor and total processor
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    //if current processor is the first one then ask for N value
    if(my_rank == 0){
        printf("Enter a N value: ");
        fflush(stdout);
        scanf("%ld", &N);
        //get current clock time
        clock_gettime(CLOCK_MONOTONIC, &start);
    }
    //broadcast N value to all other processors
    MPI_Bcast(&N, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    //calculate gap value for each processor
    localN = N / p;
    startIndex = my_rank * localN; //cal start index for each processor
    endIndex = (my_rank == p - 1)? N : startIndex + localN; //the last index needs to take all

    //for each local processor has their own start index and end index part
    for(i = startIndex; i < endIndex; i ++){
        sum += 4.0 / (1 + pow((2.0 * i + 1.0)/(2.0 * N), 2));
    }

    //combine all the local sum into total sum to get the final sum
    MPI_Reduce(&sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0 ,MPI_COMM_WORLD);
    //wait for all processes before ending
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(my_rank == 0){
        piVal = total_sum / (double)N;
        // Get the clock current time again
        // Subtract end from start to get the CPU time used.
        clock_gettime(CLOCK_MONOTONIC, &end);
        time_taken = (end.tv_sec - start.tv_sec) * 1e9;
        time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
        printf("Calculated Pi value (Parallel-AlgoI) = %12.9f\n", piVal);
        printf("Overall time (s): %lf\n", time_taken); // ts
    }

    MPI_Finalize();
    return 0;
}