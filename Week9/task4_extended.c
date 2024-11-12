#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdbool.h>
#include <time.h>

#define SHIFT_X 0
#define SHIFT_Y 1
#define SHIFT_Z 2
#define DISP 1
#define NUM_PRIMES 500

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
    int prime = 0,rec_x_lo = 0, rec_x_hi = 0, rec_y_lo = 0, rec_y_hi = 0, rec_z_lo = 0, rec_z_hi = 0;
    char filename[20];
    char *buffer, *buffer_x_lo, *buffer_x_hi, *buffer_y_lo, *buffer_y_hi, *buffer_z_lo, *buffer_z_hi; 
    int buf_size, buf_size_int, position = 0;
    int primes[NUM_PRIMES];
    int received_primes[NUM_PRIMES], received_primes_x_hi[NUM_PRIMES], received_primes_x_lo[NUM_PRIMES],
    received_primes_y_lo[NUM_PRIMES], received_primes_y_hi[NUM_PRIMES], received_primes_z_lo[NUM_PRIMES],
    received_primes_z_hi[NUM_PRIMES];

    FILE *log_file;
    struct timespec start, end, startComp, endComp;
    double time_taken;
    MPI_Request send_requests[6], recv_requests[6];

    /* start up initial MPI environment */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    clock_gettime(CLOCK_MONOTONIC, &start);

    //determine size of buffer needed to pack
    MPI_Pack_size(NUM_PRIMES, MPI_INT, MPI_COMM_WORLD, &buf_size_int);
    buf_size = buf_size_int + MPI_BSEND_OVERHEAD;
    //create buffer for all neighbours
    buffer = (char *) malloc((unsigned) buf_size);
    buffer_x_lo = (char *) malloc((unsigned) buf_size);
    buffer_x_hi = (char *) malloc((unsigned) buf_size);
    buffer_y_lo = (char *) malloc((unsigned) buf_size);
    buffer_y_hi = (char *) malloc((unsigned) buf_size);
    buffer_z_lo = (char *) malloc((unsigned) buf_size);
    buffer_z_hi = (char *) malloc((unsigned) buf_size);

    //seed the random number generator with a unique value
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
    wrap_around[0] = wrap_around[1] = wrap_around[2] =0; /* periodic shift is
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
    
    //shift the neighbours
    MPI_Cart_shift(comm3D, SHIFT_X, DISP, &nbr_x_lo, &nbr_x_hi);
    MPI_Cart_shift(comm3D, SHIFT_Y, DISP, &nbr_y_lo, &nbr_y_hi);
    MPI_Cart_shift(comm3D, SHIFT_Z, DISP, &nbr_z_lo, &nbr_z_hi);

    for(int i = 0; i < 500; i++){
        primes[i] = randomPrime(1,1000);
    }
    //pack the prime into the buffer
    MPI_Pack(primes, NUM_PRIMES, MPI_INT, buffer, buf_size, &position, MPI_COMM_WORLD);

    //send the prime buffer to neighbours
    MPI_Isend(buffer, position, MPI_PACKED, nbr_x_lo, 0, MPI_COMM_WORLD, &send_requests[0]);
    MPI_Isend(buffer, position, MPI_PACKED, nbr_x_hi, 0, MPI_COMM_WORLD, &send_requests[1]);
    MPI_Isend(buffer, position, MPI_PACKED, nbr_y_lo, 0, MPI_COMM_WORLD, &send_requests[2]);
    MPI_Isend(buffer, position, MPI_PACKED, nbr_y_hi, 0, MPI_COMM_WORLD, &send_requests[3]);
    MPI_Isend(buffer, position, MPI_PACKED, nbr_z_lo, 0, MPI_COMM_WORLD, &send_requests[4]);
    MPI_Isend(buffer, position, MPI_PACKED, nbr_z_hi, 0, MPI_COMM_WORLD, &send_requests[5]);

    //receive from neighbours with the buffer
    MPI_Irecv(buffer_x_lo, buf_size, MPI_PACKED, nbr_x_lo, 0, MPI_COMM_WORLD, &recv_requests[0]);
    MPI_Irecv(buffer_x_hi, buf_size, MPI_PACKED, nbr_x_hi, 0, MPI_COMM_WORLD, &recv_requests[1]);
    MPI_Irecv(buffer_y_lo, buf_size, MPI_PACKED, nbr_y_lo, 0, MPI_COMM_WORLD, &recv_requests[2]);
    MPI_Irecv(buffer_y_hi, buf_size, MPI_PACKED, nbr_y_hi, 0, MPI_COMM_WORLD, &recv_requests[3]);
    MPI_Irecv(buffer_y_lo, buf_size, MPI_PACKED, nbr_z_lo, 0, MPI_COMM_WORLD, &recv_requests[4]);
    MPI_Irecv(buffer_y_hi, buf_size, MPI_PACKED, nbr_z_hi, 0, MPI_COMM_WORLD, &recv_requests[5]);

    MPI_Waitall(6, send_requests, MPI_STATUSES_IGNORE);
    MPI_Waitall(6, recv_requests, MPI_STATUSES_IGNORE);

    //unpack the buffer and store in respective buffer and reset the  position
    position = 0;
    MPI_Unpack(buffer_x_lo, buf_size, &position, received_primes_x_lo, NUM_PRIMES, MPI_INT, MPI_COMM_WORLD);

    position = 0;
    MPI_Unpack(buffer_x_hi, buf_size, &position, received_primes_x_hi, NUM_PRIMES, MPI_INT, MPI_COMM_WORLD);

    position = 0;
    MPI_Unpack(buffer_y_lo, buf_size, &position, received_primes_y_lo, NUM_PRIMES, MPI_INT, MPI_COMM_WORLD);

    position = 0;
    MPI_Unpack(buffer_y_hi, buf_size, &position, received_primes_y_hi, NUM_PRIMES, MPI_INT, MPI_COMM_WORLD);

    position = 0;
    MPI_Unpack(buffer_z_lo, buf_size, &position, received_primes_z_lo, NUM_PRIMES, MPI_INT, MPI_COMM_WORLD);

    position = 0;
    MPI_Unpack(buffer_z_hi, buf_size, &position, received_primes_z_hi, NUM_PRIMES, MPI_INT, MPI_COMM_WORLD);

    sprintf(filename, "task4_ex_rank_%d.txt", my_rank);

    log_file = fopen(filename, "a");
    if (log_file == NULL){
        printf("Error opening log file\n");
        MPI_Finalize();
        return 1;
    }
    //check the prime number
    for (int i = 0; i < NUM_PRIMES; i++){
        
        if (primes[i] == received_primes_x_lo[i]) {
            fprintf(log_file, "Rank %d: Prime %d is equal to adjacent prime from rank %d.\n", my_rank, primes[i], nbr_x_lo);
        }
        if (primes[i] == received_primes_x_hi[i]) {
            fprintf(log_file, "Rank %d: Prime %d is equal to adjacent prime from rank %d.\n", my_rank, primes[i], nbr_x_hi);
        }
        if (primes[i] == received_primes_y_lo[i]) {
            fprintf(log_file, "Rank %d: Prime %d is equal to adjacent prime from rank %d.\n", my_rank, primes[i], nbr_y_lo);
        }
        if (primes[i] == received_primes_y_hi[i]) {
            fprintf(log_file, "Rank %d: Prime %d is equal to adjacent prime from rank %d.\n", my_rank, primes[i], nbr_y_hi);
        }
        if (primes[i] == received_primes_z_lo[i]) {
            fprintf(log_file, "Rank %d: Prime %d is equal to adjacent prime from rank %d.\n", my_rank, primes[i], nbr_z_lo);
        }
        if (primes[i] == received_primes_z_hi[i]) {
            fprintf(log_file, "Rank %d: Prime %d is equal to adjacent prime from rank %d.\n", my_rank, primes[i], nbr_z_hi);
        }

        printf("Prime: %d\n", primes[i]);
        printf("rank: %d, Received primes: %d, %d, %d, %d, %d, %d\n", my_cart_rank, received_primes_x_lo[i], received_primes_x_hi[i], 
                received_primes_y_lo[i], received_primes_y_hi[i], received_primes_y_lo[i], received_primes_y_hi[i]);
        

    }

        fclose(log_file);

        
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Global rank: %d. Cart rank: %d. Coord: (%d, %d).Left: %d. Right: %d. Top: %d. Bottom: %d. Front:%d. Rear:%d\n", my_rank,my_cart_rank, 
            coord[0], coord[1], nbr_x_lo, nbr_x_hi, nbr_y_lo, nbr_y_hi, nbr_z_lo, nbr_z_hi);
    
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_taken = (end.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("Overall time each processors: %f sec \n", time_taken);

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

//generate prime number
int randomPrime(int lower, int upper){
    int prime = 0;
    do {
        prime = (rand() % (upper - lower + 1)) + lower;
        } while (!isprime(prime));
    return prime;
}