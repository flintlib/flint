/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <math.h>
#include "profiler.h"
#include "fmpz.h"

#define ntests 10

void
fmpz_pow_ui_old(fmpz_t f, const fmpz_t g, ulong exp)
{
    fmpz c1;

    if (exp == WORD(0))
    {
        fmpz_one(f);
        return;
    }

    c1 = *g;

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        ulong u1 = FLINT_ABS(c1);
        ulong bits = FLINT_BIT_COUNT(u1);
        if (u1 <= UWORD(1))
        {
            fmpz_set_ui(f, u1);
        }
        else if (exp * bits <= SMALL_FMPZ_BITCOUNT_MAX)
        {
            fmpz_set_ui(f, n_pow(u1, exp));
        }
        else
        {
            __mpz_struct * mf = _fmpz_promote_val(f);

            flint_mpz_set_ui(mf, u1);
            flint_mpz_pow_ui(mf, mf, exp);
            _fmpz_demote_val(f);    /* may actually fit into a small after all */
        }

        if ((c1 < WORD(0)) && (exp & 1)) /* sign is -ve if exp odd and g -ve */
            fmpz_neg(f, f);
    }
    else
    {
        __mpz_struct * mf = _fmpz_promote_val(f);
        flint_mpz_pow_ui(mf, COEFF_TO_PTR(c1), exp);
        /* no need to demote as it can't get smaller */
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
   
    for (ix = 0; ix < 10000 * count; ix++)
    {
        fmpz_randtest(a, state, bits);
        b = n_randint(state, 40);

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_pow_ui(res, a, b);
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
   
    for (ix = 0; ix < 10000 * count; ix++)
    {
        fmpz_randtest(a, state, bits);
        b = n_randint(state, 40);

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_pow_ui_old(res, a, b);
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

    for (bits = 5; bits <= 250; bits += 5)
    {
        prof_repeat(&minnew, &maxnew, sample_new, &bits);
        prof_repeat(&minold, &maxold, sample_old, &bits);
        
        flint_printf("%d bits:      min %.2fx,    max %.2fx\n",
                bits, minold / minnew, maxold / maxnew);
    }

    return 0;
}
