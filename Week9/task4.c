#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>|
#include <stdbool.h>



#define SHIFT_X 0
#define SHIFT_Y 1
#define SHIFT_Z 2
#define DISP 1
int randomPrime(int lower, int upper);
bool isprime(int number);
int main(int argc, char *argv[]) {
    int ndims=3, size, my_rank, reorder, my_cart_rank, ierr;
    int nrows, ncols, nz;
    int nbr_x_lo, nbr_x_hi;
    int nbr_y_lo, nbr_y_hi;
    int nbr_z_lo, nbr_z_hi;
    MPI_Comm comm3D;
    int dims[ndims],coord[ndims];
    int wrap_around[ndims];
    char filename[20];
    FILE *log_file;
    int prime = 0,rec_x_lo = 0, rec_x_hi = 0, rec_y_lo = 0, rec_y_hi = 0, rec_z_lo = 0, rec_z_hi = 0;
    struct timespec start, end, startComp, endComp;
    double time_taken;


    /* start up initial MPI environment */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    srand(time(NULL) + my_rank);

    //input error checking
    if (argc == 4) {
        nrows = atoi (argv[1]);
        ncols = atoi (argv[2]);
        nz = atoi(argv[3]);
        dims[0] = nrows; /* number of rows */
        dims[1] = ncols; /* number of columns */
        dims[2] = nz;
        if( (nrows*ncols*nz) != size) {
            if( my_rank == 0) printf("ERROR: nrows*ncols*nz)=%d *%d *%d = %d != %d\n", nrows, ncols, nz, nrows*ncols*nz,size);
            MPI_Finalize();
            return 0;
        }

    } 
    
    else {
        nrows=ncols=(int)cbrt(size);
        dims[0]=dims[1]=dims[2]=0;
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
    wrap_around[0] = wrap_around[1] = wrap_around[2] = 0; /* periodic shift is
    .false. */
    reorder = 1;
    ierr =0;
    ierr= MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, wrap_around, reorder, &comm3D);
    if(ierr != 0) printf("ERROR[%d] creating CART\n",ierr);
    /* find my coordinates in the cartesian communicator group */
    MPI_Cart_coords(comm3D, my_rank, ndims, coord);
    /* use my cartesian coordinates to find my rank in cartesian
    group*/
    MPI_Cart_rank(comm3D, coord, &my_cart_rank);
    /* get my neighbors; axis is coordinate dimension of shift */
    /* axis=0 ==> shift along the rows: P[my_row-1]: P[me] :
    P[my_row+1] */
    /* axis=1 ==> shift along the columns P[my_col-1]: P[me] :
    P[my_col+1] */
    
    //shift the cart to let the processor know the neighbours
    MPI_Cart_shift(comm3D, SHIFT_X, DISP, &nbr_x_lo, &nbr_x_hi);
    MPI_Cart_shift(comm3D, SHIFT_Y, DISP, &nbr_y_lo, &nbr_y_hi);
    MPI_Cart_shift(comm3D, SHIFT_Z, DISP, &nbr_z_lo, &nbr_z_hi);

    //iterate 500 times
    for(int i = 0; i < 500; i++){
    
        prime = randomPrime(1, 1000); //get a random prime number
    
        //mpi send to its neighbours
        MPI_Send(&prime, 1, MPI_INT, nbr_x_lo, 0, MPI_COMM_WORLD);
        MPI_Send(&prime, 1, MPI_INT, nbr_x_hi, 0, MPI_COMM_WORLD);
        MPI_Send(&prime, 1, MPI_INT, nbr_y_lo, 0, MPI_COMM_WORLD);
        MPI_Send(&prime, 1, MPI_INT, nbr_y_hi, 0, MPI_COMM_WORLD);
        MPI_Send(&prime, 1, MPI_INT, nbr_z_lo, 0, MPI_COMM_WORLD);
        MPI_Send(&prime, 1, MPI_INT, nbr_z_hi, 0, MPI_COMM_WORLD);

        //receive from its neighbour
        MPI_Recv(&rec_x_lo, 1, MPI_INT, nbr_x_lo, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&rec_x_hi, 1, MPI_INT, nbr_x_hi, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&rec_y_lo, 1, MPI_INT, nbr_y_lo, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&rec_y_hi, 1, MPI_INT, nbr_y_hi, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&rec_z_hi, 1, MPI_INT, nbr_z_lo, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&rec_z_hi, 1, MPI_INT, nbr_z_hi, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //determine the file name
        sprintf(filename, "task4_rank_%d.txt", my_rank);
        //open the filename
        log_file = fopen(filename, "a");
        if (log_file == NULL){
            printf("Error opening log file\n");
            MPI_Finalize();
            return 1;
        }
        //if it matches the prime then write in file
        if (prime == rec_x_lo) {
            fprintf(log_file, "Rank %d: Prime %d is equal to adjacent prime from rank %d.\n", my_rank, prime, nbr_x_lo);
        }
        if (prime == rec_x_hi) {
            fprintf(log_file, "Rank %d: Prime %d is equal to adjacent prime from rank %d.\n", my_rank, prime, nbr_x_hi);
        }
        if (prime == rec_y_lo) {
            fprintf(log_file, "Rank %d: Prime %d is equal to adjacent prime from rank %d.\n", my_rank, prime, nbr_y_lo);
        }
        if (prime == rec_y_hi) {
            fprintf(log_file, "Rank %d: Prime %d is equal to adjacent prime from rank %d.\n", my_rank, prime, nbr_y_hi);
        }
        if (prime == rec_z_lo) {
            fprintf(log_file, "Rank %d: Prime %d is equal to adjacent prime from rank %d.\n", my_rank, prime, nbr_z_lo);
        }
        if (prime == rec_z_hi) {
            fprintf(log_file, "Rank %d: Prime %d is equal to adjacent prime from rank %d.\n", my_rank, prime, nbr_z_hi);
        }

        fclose(log_file);

        printf("Prime: %d\n", prime);

        printf("rank: %d, Received primes: %d, %d, %d, %d, %d, %d\n", 
        my_cart_rank, rec_x_lo, rec_x_hi, rec_y_lo, rec_y_hi, rec_z_lo, rec_z_hi);
    }
    //wait for all processes
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Global rank: %d, Cartesian rank: %d, Coordinates: (%d, %d, %d)\n",
           my_rank, my_cart_rank, coord[0], coord[1], coord[2]);

    printf("Neighbors: Left: %d, Right: %d, Top: %d, Bottom: %d, Front: %d, Rear: %d\n",
           nbr_x_lo, nbr_x_hi, nbr_y_lo, nbr_y_hi, nbr_z_lo, nbr_z_hi);
           
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_taken = (end.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("Overall time each processors: %f sec \n", time_taken);

    fflush(stdout);
    MPI_Comm_free( &comm3D );
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

//generate random prime number
int randomPrime(int lower, int upper){
    int prime = 0;
    do {
        prime = (rand() % (upper - lower + 1)) + lower;
        } while (!isprime(prime));
    return prime;
}