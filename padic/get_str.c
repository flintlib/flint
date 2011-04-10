#include <stdlib.h>
#include <stdio.h>

#include <limits.h>

#include "padic.h"

/*
    Returns the number of digits in the base $b$ representation 
    of the integer $n$.

    Assumes that $b \geq 2$.
 */
static long sizeinbase_si(long n, long b)
{
    long c;

    if (n > 0)
        c = 0;
    else if (n > LONG_MIN)
    {
        n = -n;
        c = 1;
    }
    else  /* n == LONG_MIN */
    {
        if (n % b)
        {
            n = - (n + 1);
            c = 1;
        }
        else
        {
            n = - (n / b);
            c = 2;
        }
    }

    for ( ; n > 0; n /= b, c++) ;

    return c;
}

char * padic_get_str(const padic_t op, const padic_ctx_t ctx)
{
    return NULL;
}

