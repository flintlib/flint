/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "profiler.h"
#include "fmpz.h"

#define ntests 1000

void
fmpz_mul_2exp_old(fmpz_t f, const fmpz_t g, ulong exp)
{
    fmpz d = *g;

    if (!COEFF_IS_MPZ(d))       /* g is small */
    {
        ulong dabs = FLINT_ABS(d);
        ulong bits = FLINT_BIT_COUNT(dabs);
        if (bits == 0)
        {
            fmpz_set_si(f, 0);
        }
        else if (bits + exp <= SMALL_FMPZ_BITCOUNT_MAX)  /* result will fit in small */
        {
            fmpz_set_si(f, d << exp);
        }
        else                    /* result is large */
        {
            __mpz_struct * mf = _fmpz_promote(f);   /* g is saved */
            flint_mpz_set_si(mf, d);
            mpz_mul_2exp(mf, mf, exp);
        }
    }
    else                        /* g is large */
    {
        __mpz_struct * mf = _fmpz_promote(f);   /* g is already large */
        mpz_mul_2exp(mf, COEFF_TO_PTR(d), exp);
    }
}

void
sample_new(void * arg, ulong count)
{
    fmpz_t res, a;
    ulong ix, jx, b;
    int bits = *((int *) arg);

    FLINT_TEST_INIT(state);

    fmpz_init(res);
    fmpz_init(a);
   
    for (ix = 0; ix < count; ix++)
    {
        fmpz_randtest(a, state, bits);
        b = n_randint(state, 200);

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_mul_2exp(res, a, b);
        prof_stop();
    }

    flint_randclear(state);
    fmpz_clear(res);
    fmpz_clear(a);
}

void
sample_old(void * arg, ulong count)
{
    fmpz_t res, a;
    ulong ix, jx, b;
    int bits = *((int *) arg);

    FLINT_TEST_INIT(state);

    fmpz_init(res);
    fmpz_init(a);
   
    for (ix = 0; ix < count; ix++)
    {
        fmpz_randtest(a, state, bits);
        b = n_randint(state, 200);

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_mul_2exp_old(res, a, b);
        prof_stop();
    }

    flint_randclear(state);
    fmpz_clear(res);
    fmpz_clear(a);
}

int
main(void)
{
    double minnew, maxnew, minold, maxold;
    int bits;

    for (bits = 5; bits <= 150; bits += 5)
    {
        prof_repeat(&minnew, &maxnew, sample_new, &bits);
        prof_repeat(&minold, &maxold, sample_old, &bits);
        
        flint_printf("%d bits:      min %.2fx,    max %.2fx\n",
                bits, minold / minnew, maxold / maxnew);
    }

    return 0;
}
