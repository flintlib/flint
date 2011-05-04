#include <stdlib.h>
#include <stdio.h>

#include "fmpz_poly_q.h"

int fmpz_poly_q_print(const fmpz_poly_q_t op)
{
    char *str;

    str = fmpz_poly_q_get_str(op);
    printf("%s", str);
    free(str);

    return 1;
}
