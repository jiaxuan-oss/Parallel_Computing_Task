#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
int main(int argc, char* argv[])
{
    int my_rank;
    int p;

    //initialise MPI
    MPI_Init(&argc, &argv);
    //get current processor
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    //get all processor
    MPI_Comm_size(MPI_COMM_WORLD, &p); 
    int val = -1;

    do{
    if(my_rank == 0)
    {
    printf("Enter a round number (> 0 ): ");
    fflush(stdout);
    scanf("%d", &val);
    }

    //broadcast to all the processors
    // Please insert one line of code here
    MPI_Bcast(&val, 1, MPI_INT, 0, MPI_COMM_WORLD);
    printf("Processors: %d. Received Value: %d\n", my_rank, val);
    fflush(stdout);
    }while(val > 0);

    MPI_Finalize();
    return 0;
}