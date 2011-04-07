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

int padic_fprint(FILE * file, const padic_t op, const padic_ctx_t ctx)
{
    if (padic_is_zero(op, ctx))
    {
        fprintf(file, "0 + O(");
        fmpz_fprint(file, ctx->p);
        fprintf(file, "^%ld)", ctx->N);

        return 1;
    }
    
    if (ctx->mode == PADIC_TERSE)
    {
        if (op[1] >= 0)
        {
            fmpz_t pow, x;

            fmpz_init(pow);
            fmpz_init(x);

            fmpz_pow_ui(pow, ctx->p, op[1]);
            fmpz_mul(x, pow, op);

            fmpz_fprint(file, x);

            fmpz_clear(pow);
            fmpz_clear(x);
        }
        else
        {
            fmpz_t pow;
            mpq_t y;

            fmpz_init(pow);
            mpq_init(y);

            fmpz_pow_ui(pow, ctx->p, op[1]);
            fmpz_get_mpz(mpq_numref(y), op);
            fmpz_get_mpz(mpq_denref(y), pow);

            gmp_fprintf(file, "%Qd", y);

            fmpz_clear(pow);
            mpq_clear(y);
        }
    }
    else if (ctx->mode == PADIC_SERIES)
    {
        fmpz_t x;
        fmpz_t d;
        long j;

        fmpz_init(d);
        fmpz_init(x);

        fmpz_set(x, op);
        
        for (j = 0; j < ctx->N - op[1]; j++)
        {
            fmpz_mod(d, x, ctx->p);       /* d = u mod p^{j+1} */
            fmpz_sub(x, x, d);            /* x = x - d */
            fmpz_divexact(x, x, ctx->p);  /* x = x / p */

            if (!fmpz_is_zero(d))
            {
                if(j + op[1] != 0)
                {
                    fmpz_fprint(file, d);
                    fprintf(file, "*");
                    fmpz_fprint(file, ctx->p);
                    fprintf(file, "^%ld + ", j + op[1]);
                }
                else
                {
                    fmpz_fprint(file, d);
                    fprintf(file, " + ");
                }
            }
        }

        fprintf(file, "O(");
        fmpz_fprint(file, ctx->p);
        fprintf(file, "^%ld)", ctx->N);

        fmpz_clear(x);
        fmpz_clear(d);
    }
    else if (ctx->mode == PADIC_VAL_UNIT)
    {
        fmpz_fprint(file, op);
        fprintf(file, "*");
        fmpz_fprint(file, ctx->p);
        fprintf(file, "^%ld + O(", op[1]);
        fmpz_fprint(file, ctx->p);
        fprintf(file, "^%ld)", ctx->N);
    }
    else
    {
        printf("Exception (padic_fprint).  Unknown printing mode.\n");
        abort();
    }

    return 1;
}
