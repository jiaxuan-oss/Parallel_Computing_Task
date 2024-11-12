#include <stdio.h>
#include <math.h>

int main() {
	int counter = 0;
	printf("Using for loop");
	//loop 1: For loop
	//for loop print until the condition <300 is met
	for(int i = 0; i < 300; i = i + 1){
		printf("Using for loop\n");

	}
	//Loop 2: while loop
	//while loop print until the counter is <300
	printf("Using While loop");
	while (counter < 300){
		printf("Using While loop\n");
		counter = counter + 1;
	
	}
	return (0);
}
