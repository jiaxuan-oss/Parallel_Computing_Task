#include <stdio.h>
void swap(int *a, int*b);
void printArray(int *arr, int size);

void main(){
    //initialise and declare the list
    int data[] = {1,2,3,4,5};;
    //call the function of printarray
    //function 2 with pointer
    printArray(data, 5);

}


void swap(int *a, int *b){
    //declaration
    int temp;
    //swapping value with temp variable
     
    temp = *a; //temp now holding value of address a
    *a = *b;   //replacing value of address a with value of address b 
    *b = temp; //value of addresss b is replaced with temp
}

void printArray(int *arr, int size){
    //looping the array elem
    for (int i = 0; i < size; i++){
        //if is not the last elem 
        if(i < size - 1){
            //then swap
            //Function 1 with pointer
            swap(&arr[i], &arr[i+1]);
        }
    }
    //using for loop to print the value in the array
    for (int i = 0; i < size; i++){
        printf("%d\n",arr[i]);    
        
    }

}

