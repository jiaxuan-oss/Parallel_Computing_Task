#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <pthread.h>
#define THREADS 8

bool isprime(int number);

//define global variable
size_t num_elem = 10000000;
int *arr;
int counter = 0;
pthread_mutex_t mutex;

//call prime number function 
void *calprime (void *arg){

    //indicate start and end
    int start = *((int*)arg);
    int end = *((int*)arg + 1); 

    //looping through the numbers to check prime
    for (int i = start; i < end ; i++){
        if (isprime(i)){

            //if is prime then lock the threads so that one thread can update each time
            //to prevent race condition
            pthread_mutex_lock(&mutex);
            arr[counter] = i;
            counter++;
            pthread_mutex_unlock(&mutex);
            //unlock after done updating
        }
    }
}

int main () 
{
    //declaration and initialisation
    pthread_mutex_init(&mutex, NULL);
    struct timespec start, end, startComp, endComp;
    double time_taken;
    pthread_t threads[THREADS];
    int thread_args[THREADS][2];
    pthread_mutex_t mutex;
    arr = (int *)malloc(num_elem * sizeof(int));

    //get step for each threads
    int step = 10000000 / THREADS;

     // Get current clock time.
    printf("Compute:\n");
    clock_gettime(CLOCK_MONOTONIC, &start);

    for(int i = 0; i < THREADS; i++){
        //set each thread for its start and end in a 2d array
        thread_args[i][0] = i * step;
        if (i == THREADS - 1){
            thread_args[i][1] = 10000000;
        }
        else{
            thread_args[i][1] = (i + 1) * step;
        }
        //create the threads using i 
        pthread_create(&threads[i], NULL, calprime, &thread_args[i]);
    }

    //join the threads together
    for (int i = 0; i < THREADS; i++){
        pthread_join(threads[i], NULL);
    }

    pthread_mutex_destroy(&mutex);
    clock_gettime(CLOCK_MONOTONIC, &end);
    //print out the updated array-> all the prime number
    for(int i = 0; i < counter; i++){
        printf("%d ",arr[i]);
    }

    
    //get the clock time and get the time that the application run
    time_taken = (end.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("Overall time: %f sec \n", time_taken);

}

//determine prime number
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