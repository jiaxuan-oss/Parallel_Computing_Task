#include <stdio.h>
#include <math.h>

int main(){
    //declaration and initialisation
    int flag = 1;
    int counter = 0;
    int day = 2;
    //while loop until the flag is 0
    while (flag){
		counter = counter + 1;
        //conditional 1: if else, if the counter == 10 then stop the loop
        if(counter == 10){
            flag = 0;
        }

        else{
            printf("Havent ready yet\n");
        }
	}
    //conditional 2: Switch case
    //using switch case to determine day
    //if is 1 then monday 2 tuesday 3 wednesday...
    switch (day) {
        case 1:
            printf("Monday\n");
            break;
        case 2:
            printf("Tuesday\n");
            break;
        case 3:
            printf("Wednesday\n");
            break;
        case 4:
            printf("Thursday\n");
            break;
        case 5:
            printf("Friday\n");
            break;
        case 6:
            printf("Saturday\n");
            break;
        case 7:
            printf("Sunday\n");
            break;
        default:
            printf("Invalid day\n");
            break;
    }

    return (0);
}