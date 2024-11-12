#include <stdio.h>
#include <ctype.h>

int main()
{
    //declare and initialisation 
    int c, word_count = 0;
    int current_word = 0;

    while ( (c = getchar()) != EOF )
    {   
        //check c is a newline
        if (c == '\n') 
        {
            current_word = 0;
        }
        //check c is a space
        else if (isspace(c)) 
        { 
            //if is a space then current word = 0
            current_word = 0;   
        }
        else 
        {
            if (!current_word)
            {   
                //if is a character then -> a new word and current word is 1
                current_word = 1;
                word_count++; //increment word count
            }

        }
    }

    printf("number of words = %d\n", word_count);

    return(0);
}
