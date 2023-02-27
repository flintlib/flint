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

#define ntests 2000

void
fmpz_mul_ui_old(fmpz_t f, const fmpz_t g, ulong x)
{
    fmpz c2 = *g;

    if (!COEFF_IS_MPZ(c2)) /* c2 is small */
    {
        mp_limb_t th, tl;
        mp_limb_t uc2 = FLINT_ABS(c2);

        /* unsigned limb by limb multiply (assembly for most CPU's) */
        umul_ppmm(th, tl, uc2, x);
        if (c2 >= 0)
            fmpz_set_uiui(f, th, tl);
        else
            fmpz_neg_uiui(f, th, tl);
    }
    else                        /* c2 is large */
    {
        if (x == 0)
            fmpz_zero(f);
        else
        {
            /* Promote without val as if aliased both are large */
            __mpz_struct *mpz_ptr = _fmpz_promote(f);
            mpz_mul_ui(mpz_ptr, COEFF_TO_PTR(c2), x);
        }
    }
}

void
sample_new(void * arg, ulong count)
{
    fmpz_t res, a;
    ulong b;
    ulong ix, jx;
    int bits = *((int *) arg);

    FLINT_TEST_INIT(state);

    fmpz_init(res);
    fmpz_init(a);
   
    for (ix = 0; ix < count; ix++)
    {
        fmpz_randtest(a, state, bits + 20);
        b = n_randbits(state, bits);

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_mul_ui(res, a, b);
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
    ulong b;
    ulong ix, jx;
    int bits = *((int *) arg);

    FLINT_TEST_INIT(state);

    fmpz_init(res);
    fmpz_init(a);
   
    for (ix = 0; ix < count; ix++)
    {
        fmpz_randtest(a, state, bits + 20);
        b = n_randbits(state, bits);

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_mul_ui_old(res, a, b);
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

    for (bits = 1; bits <= 64; bits += 3)
    {
        prof_repeat(&minnew, &maxnew, sample_new, &bits);
        prof_repeat(&minold, &maxold, sample_old, &bits);
        
        flint_printf("%d bits:      min %.2fx,    max %.2fx\n",
                bits, minold / minnew, maxold / maxnew);
    }

    return 0;
}
