/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "gmpcompat.h"

#define ntests 30

void
fmpz_mul_old(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1, c2;
    __mpz_struct *z;

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

    z = _fmpz_promote(f); /* h is saved, g is already large */

    if (!COEFF_IS_MPZ(c2))      /* g is large, h is small */
        flint_mpz_mul_si(z, COEFF_TO_PTR(c1), c2);
    else                        /* c1 and c2 are large */
        mpz_mul(z, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
}

void
sample_new(void * arg, ulong count)
{
    fmpz *res, *a, *b;
    ulong ix, jx;
    int bits = *((int *) arg);

    FLINT_TEST_INIT(state);

    res = _fmpz_vec_init(ntests);
    a = _fmpz_vec_init(ntests);
    b = _fmpz_vec_init(ntests);

    for (ix = 0; ix < 10 * count; ix++)
    {
        for (jx = 0; jx < ntests; jx++)
        {
            fmpz_randtest(a + jx, state, bits);
            fmpz_randtest(b + jx, state, bits);
            fmpz_randtest(res + jx, state, bits);
        }

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_mul(res + jx, a + jx, b + jx);
        prof_stop();
    }

    _fmpz_vec_clear(res, ntests);
    _fmpz_vec_clear(a, ntests);
    _fmpz_vec_clear(b, ntests);
    flint_randclear(state);
}

void
sample_old(void * arg, ulong count)
{
    fmpz *res, *a, *b;
    ulong ix, jx;
    int bits = *((int *) arg);

    FLINT_TEST_INIT(state);

    res = _fmpz_vec_init(ntests);
    a = _fmpz_vec_init(ntests);
    b = _fmpz_vec_init(ntests);

    for (ix = 0; ix < 10 * count; ix++)
    {
        for (jx = 0; jx < ntests; jx++)
        {
            fmpz_randtest(a + jx, state, bits);
            fmpz_randtest(b + jx, state, bits);
            fmpz_randtest(res + jx, state, bits);
        }

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_mul_old(res + jx, a + jx, b + jx);
        prof_stop();
    }

    _fmpz_vec_clear(res, ntests);
    _fmpz_vec_clear(a, ntests);
    _fmpz_vec_clear(b, ntests);
    flint_randclear(state);
}

slong sizes[] = { 10, 30, 60, 62, 64, 66, 80, 128, 160, 256, 512, 1024, 4096, 0 };

int
main(void)
{
    double minnew, maxnew, minold, maxold;
    int i, bits;

    for (i = 0; (bits = sizes[i]) != 0; i++)
    {
        prof_repeat(&minnew, &maxnew, sample_new, &bits);
        prof_repeat(&minold, &maxold, sample_old, &bits);

        flint_printf("%d bits:      min %.2fx,    max %.2fx\n",
                bits, minold / minnew, maxold / maxnew);
    }

    return 0;
}
