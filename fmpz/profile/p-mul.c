/*
    Copyright (C) 2021 Albin Ahlbäck

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
fmpz_mul_old(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1, c2;
    __mpz_struct *mpz_ptr;

    c1 = *g;

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        fmpz_mul_si(f, h, c1);
        return;
    }

    c2 = *h;                    /* save h in case it is aliased with f */

    if (c2 == WORD(0))               /* special case, h = 0  */
    {
        fmpz_zero(f);
        return;
    }

    mpz_ptr = _fmpz_promote(f); /* h is saved, g is already large */

    if (!COEFF_IS_MPZ(c2))      /* g is large, h is small */
        flint_mpz_mul_si(mpz_ptr, COEFF_TO_PTR(c1), c2);
    else                        /* c1 and c2 are large */
        mpz_mul(mpz_ptr, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
}

void
sample_new(void * arg, ulong count)
{
    fmpz_t res, a, b;
    ulong ix, jx;
    int bits = *((int *) arg);

    FLINT_TEST_INIT(state);

    fmpz_init(res);
    fmpz_init(a);
    fmpz_init(b);
   
    for (ix = 0; ix < count; ix++)
    {
        fmpz_randtest(a, state, bits);
        fmpz_randtest(b, state, bits);

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_mul(res, a, b);
        prof_stop();
    }

    flint_randclear(state);
    fmpz_clear(res);
    fmpz_clear(a);
    fmpz_clear(b);
}

void
sample_old(void * arg, ulong count)
{
    fmpz_t res, a, b;
    ulong ix, jx;
    int bits = *((int *) arg);

    FLINT_TEST_INIT(state);

    fmpz_init(res);
    fmpz_init(a);
    fmpz_init(b);
   
    for (ix = 0; ix < count; ix++)
    {
        fmpz_randtest(a, state, bits);
        fmpz_randtest(b, state, bits);

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_mul_old(res, a, b);
        prof_stop();
    }

    flint_randclear(state);
    fmpz_clear(res);
    fmpz_clear(a);
    fmpz_clear(b);
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
