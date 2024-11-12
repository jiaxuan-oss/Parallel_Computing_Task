#include <stdio.h>
#include <mpi.h>

//structure definition
struct valuestruct {
    int a;
    double b;
};

int main(int argc, char** argv)
{
    //variable declarations
    struct valuestruct values;
    int myrank;
    MPI_Datatype Valuetype; //custom mpi datatype
    MPI_Datatype type[2] = { MPI_INT, MPI_DOUBLE }; 
    int blocklen[2] = { 1, 1};
    MPI_Aint disp[2]; //array to store offset of each member

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    //get memory addresses of value.a and b and store in disp
    MPI_Get_address(&values.a, &disp[0]);
    MPI_Get_address(&values.b, &disp[1]);
    //Make relative, reflect the offsets within the structure
    //know how many step to jump from first address to second address
    //so store the displacement
    disp[1]=disp[1]-disp[0];
    disp[0]=0;

    // Create MPI struct
    // Please insert one line of code here (Hint: MPI_Type_create_struct)
    //create custom MPI datatype
    //element in struct, specifiying number of elem of each type, disp of the address, specify MPI datatype of each block, MPI custom type
    MPI_Type_create_struct(2, blocklen, disp, type, &Valuetype);
    
    //commit the custom datatype
    MPI_Type_commit(&Valuetype);

    do{
    if (myrank == 0){
    printf("Enter an round number (>0) & a real number: ");
    fflush(stdout);
    scanf("%d%lf", &values.a, &values.b);
    }
    // Please insert one line of code here (Hint: MPI_Bcast)
    //broadcast the value array to other processor
    MPI_Bcast(&values, 1, Valuetype, 0, MPI_COMM_WORLD);

    printf("Rank: %d. values.a = %d. values.b = %lf\n",
    myrank, values.a, values.b);
    fflush(stdout);
    }while(values.a > 0);
    /* Clean up the type */
    MPI_Type_free(&Valuetype);
    MPI_Finalize();

    return 0;
}