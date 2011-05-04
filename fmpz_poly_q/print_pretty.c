#include <stdlib.h>
#include <stdio.h>

#include "fmpz_poly_q.h"

int fmpz_poly_q_print_pretty(const fmpz_poly_q_t op, const char *x)
{
    char *str;

    str = fmpz_poly_q_get_str_pretty(op, x);
    printf("%s", str);
    free(str);

    return 1;
}
