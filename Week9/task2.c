#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdbool.h>
#include <time.h>

#define SHIFT_ROW 0
#define SHIFT_COL 1
#define DISP 1
int randomNum(int lower, int upper);
bool isprime(int number);
int main(int argc, char *argv[]) {
    int ndims=2, size, my_rank, reorder, my_cart_rank, ierr;
    int nrows, ncols;
    int nbr_i_lo, nbr_i_hi;
    int nbr_j_lo, nbr_j_hi;
    MPI_Comm comm2D;
    int dims[ndims],coord[ndims];
    int wrap_around[ndims];
    int prime = 0;
    int rec_i_lo = 0, rec_i_hi = 0, rec_j_lo = 0, rec_j_hi = 0;
    char filename[20];
    FILE *log_file;
    struct timespec start, end, startComp, endComp;
    double time_taken;

    /* start up initial MPI environment */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    clock_gettime(CLOCK_MONOTONIC, &start);

    //seed the random number generator with a unique value
    srand(time(NULL) + my_rank);

    /* process command line arguments*/
    if (argc == 3) {
        nrows = atoi (argv[1]);
        ncols = atoi (argv[2]);
        dims[0] = nrows; /* number of rows */
        dims[1] = ncols; /* number of columns */
        if( (nrows*ncols) != size) {
            if( my_rank == 0) printf("ERROR: nrows*ncols)=%d *%d = %d != %d\n", nrows, ncols, nrows*ncols,size);
            MPI_Finalize();
            return 0;
        }

    } 
    
    else {
        nrows=ncols=(int)sqrt(size);
        dims[0]=dims[1]=0;
    }
    /************************************************************
    */
    /* create cartesian topology for processes */
    /************************************************************
    */
    MPI_Dims_create(size, ndims, dims);
    if(my_rank == 0){
        printf("Root Rank: %d. Comm Size: %d: Grid Dimension =[%d x %d] \n",my_rank,size,dims[0],dims[1]);
    }
    
    /* create cartesian mapping */
    wrap_around[0] = wrap_around[1] = 0; /* periodic shift is
    .false. */
    reorder = 1;
    ierr =0;
    ierr= MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, wrap_around, reorder, &comm2D);
    if(ierr != 0) printf("ERROR[%d] creating CART\n",ierr);
    /* find my coordinates in the cartesian communicator group */
    MPI_Cart_coords(comm2D, my_rank, ndims, coord);
    /* use my cartesian coordinates to find my rank in cartesian
    group*/
    MPI_Cart_rank(comm2D, coord, &my_cart_rank);
    /* get my neighbors; axis is coordinate dimension of shift */
    /* axis=0 ==> shift along the rows: P[my_row-1]: P[me] :
    P[my_row+1] */
    /* axis=1 ==> shift along the columns P[my_col-1]: P[me] :
    P[my_col+1] */
    
    MPI_Cart_shift(comm2D, SHIFT_ROW, DISP, &nbr_i_lo, &nbr_i_hi);
    MPI_Cart_shift(comm2D, SHIFT_COL, DISP, &nbr_j_lo, &nbr_j_hi);

    //iterate 500
    for(int i = 0; i < 500; i++){
        do {
            prime = randomNum(1, 1000); //get prime number
        } while (!isprime(prime));
        //send prime to its neighbour
        MPI_Send(&prime, 1, MPI_INT, nbr_i_lo, 0, MPI_COMM_WORLD);
        MPI_Send(&prime, 1, MPI_INT, nbr_i_hi, 0, MPI_COMM_WORLD);
        MPI_Send(&prime, 1, MPI_INT, nbr_j_lo, 0, MPI_COMM_WORLD);
        MPI_Send(&prime, 1, MPI_INT, nbr_j_hi, 0, MPI_COMM_WORLD);

        //receive from its neighbour and store to rec 
        MPI_Recv(&rec_i_lo, 1, MPI_INT, nbr_i_lo, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&rec_i_hi, 1, MPI_INT, nbr_i_hi, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&rec_j_lo, 1, MPI_INT, nbr_j_lo, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&rec_j_hi, 1, MPI_INT, nbr_j_hi, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        sprintf(filename, "task2_rank_%d.txt", my_rank);

        log_file = fopen(filename, "a");
        if (log_file == NULL){
            printf("Error opening log file\n");
            MPI_Finalize();
            return 1;
        }
        //compare prime with neighbour and write down if is same
        if (prime == rec_i_lo) {
            fprintf(log_file, "Rank %d: Prime %d is equal to adjacent prime from rank %d.\n", my_rank, prime, nbr_i_lo);
        }
        if (prime == rec_i_hi) {
            fprintf(log_file, "Rank %d: Prime %d is equal to adjacent prime from rank %d.\n", my_rank, prime, nbr_i_hi);
        }
        if (prime == rec_j_lo) {
            fprintf(log_file, "Rank %d: Prime %d is equal to adjacent prime from rank %d.\n", my_rank, prime, nbr_j_lo);
        }
        if (prime == rec_j_hi) {
            fprintf(log_file, "Rank %d: Prime %d is equal to adjacent prime from rank %d.\n", my_rank, prime, nbr_j_hi);
        }

        fclose(log_file);

        printf("Prime: %d\n", prime);
        printf("rank: %d, Received primes: %d, %d, %d, %d\n", my_cart_rank, rec_i_lo, rec_i_hi, rec_j_lo, rec_j_hi);

    }
    //wait for all processes
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Global rank: %d. Cart rank: %d. Coord: (%d, %d).Left: %d. Right: %d. Top: %d. Bottom: %d\n", my_rank,my_cart_rank, coord[0], coord[1], nbr_j_lo, nbr_j_hi, nbr_i_lo, nbr_i_hi);
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_taken = (end.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("Overall time each processors: %f sec \n", time_taken);
    fflush(stdout);
    
    MPI_Comm_free( &comm2D );
    MPI_Finalize();
    return 0;
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

//generate random number
int randomNum(int lower, int upper){
    return (rand() % (upper - lower + 1)) + lower;
}