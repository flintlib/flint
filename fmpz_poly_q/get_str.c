#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "fmpz_poly_q.h"

/**
 * \ingroup  StringConversions
 * 
 * Returns the string representation of the rational function \c op.
 */
char * fmpz_poly_q_get_str(const fmpz_poly_q_t op)
{
    int i, j;
    char * str;
    char * numstr;
    char * denstr;
    
    if (fmpz_poly_is_one(op->den))
    {
        numstr = fmpz_poly_get_str(op->num);
        i = strlen(numstr) - 1;
        if (numstr[i] == ' ')
        {
            numstr[i] = '\0';
        }
        return numstr;
    }
    
    numstr = fmpz_poly_get_str(op->num);
    denstr = fmpz_poly_get_str(op->den);
    
    i = strlen(numstr) - 1;
    if (numstr[i] == ' ')
        numstr[i] = '\0';
    i = strlen(denstr) - 1;
    if (denstr[i] == ' ')
        denstr[i] = '\0';
    
    str = malloc(strlen(numstr) + strlen(denstr) + 2);
    if (str == NULL)
    {
        printf("ERROR (fmpz_poly_q_get_str).  Memory allocation failed.\n");
        abort();
    }
    
    for (i = 0; i < strlen(numstr); i++)
        str[i] = numstr[i];
    str[i++] = '/';
    for (j = 0; j < strlen(denstr); j++)
        str[i++] = denstr[j];
    str[i] = '\0';
    
    free(numstr);
    free(denstr);
    
    return str;
}
