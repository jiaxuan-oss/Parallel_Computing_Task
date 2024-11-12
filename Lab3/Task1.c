#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

bool isprime(int number);
int main () 
{   
    //declaration and initialisation
    int *arr = (int *)malloc(1000000 * sizeof(int));
    struct timespec start, end, startComp, endComp;
    double time_taken;
    int n;
    int counter = 0;

    // Get current clock time.

    printf("Compute:\n");
    clock_gettime(CLOCK_MONOTONIC, &start);

    //loop through the numbers to check prime
    for (int i = 0; i < 1000000 ; i++){
        if (isprime(i)){
            //update the array if there is prime
            counter++;
            arr[counter] = i;
        }
    }


    //print out the numbers in array
    for(int i = 0; i <= counter; i++){
        printf("%d\n", arr[i]);
    }
    //get time to calculate how long the algo runs
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

    //calculate prime number return false if it can be divided by a number before it
    for (int i = 3; i <= sqrt(number); i++ ){
        if (number % i == 0){
            return false;
        }
    }
    return true;

}