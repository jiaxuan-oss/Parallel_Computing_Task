#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>

bool isprime(int number);

int main () 
{   
    int *arr = (int *)malloc(1000000 * sizeof(int));
    struct timespec start, end, startComp, endComp;
    double time_taken;
    int n;
    int counter = 0;

    // Get current clock time.
    printf("Compute:\n");
    clock_gettime(CLOCK_MONOTONIC, &start);
    
    //parallel this for loop
    #pragma omp parallel
    {   
        //define local array so each thread will create these
        int *localarr = (int *)malloc(1000000 * sizeof(int));
        int localcounter = 0;

        #pragma omp for schedule(guided)
        {   
            //each thread will update their prime number in their own local array 
            for (int i = 2; i < 1000000 ; i++){
                if (isprime(i)){
                    localarr[localcounter] = i;  
                    localcounter++;
                }  
            }
        }

        //use to lock the shared memory to prevent race condition and data corruption
        #pragma omp critical
        {   
            //after they calculate they own value in their local array they update this global array accordingly
            //only one thread able to update everytime
            for(int i = 0; i < localcounter; i++){
                arr[counter] = localarr[i];
                counter++;
            }
        }
    }
    
    //print all the output
    for(int i = 0; i < counter; i++){
        printf("%d\n", arr[i]);
    }

    //get the time 
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_taken = (end.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("Overall time: %f sec \n", time_taken);

}



bool isprime(int number){

    if (number <= 1){
        return false;
    }

    if (number == 2){
        return true;
    }

    if (number % 2 == 0){
        return false;
    }
    
    //loop until sqrt of n as p*q = n so loop until one of them is sufficient
    for (int i = 3; i <= sqrt(number); i++ ){
        if (number % i == 0){
            return false;
        }
    }
    return true;

}