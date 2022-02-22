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

#define ntests 10

void fmpz_add_ui_old(fmpz_t f, const fmpz_t g, ulong x)
{
    fmpz c = *g;

    if (!COEFF_IS_MPZ(c))  /* g is small */
    {
        mp_limb_t sum[2];
        if (c >= WORD(0))  /* both operands non-negative */
        {
            add_ssaaaa(sum[1], sum[0], 0, c, 0, x);
            fmpz_set_uiui(f, sum[1], sum[0]);
        }
        else  /* coeff is negative, x positive */
        {
            if (-c > x)
                fmpz_set_si(f, x + c); /* can't overflow as g is small and x smaller */
            else
                fmpz_set_ui(f, x + c);  /* won't be negative and has to be less than x */
        }
    }
    else
    {	
        __mpz_struct * mf = _fmpz_promote(f);  /* g is already large */
        __mpz_struct * mc = COEFF_TO_PTR(c);
        flint_mpz_add_ui(mf, mc, x);
        _fmpz_demote_val(f);  /* cancellation may have occurred */
    }
}

void
fmpz_sub_ui_old(fmpz_t f, const fmpz_t g, ulong x)
{
    fmpz c = *g;

    if (!COEFF_IS_MPZ(c))       /* coeff is small */
    {
        mp_limb_t sum[2];
        if (c < WORD(0))             /* g negative, x positive, so difference is negative */
        {
            add_ssaaaa(sum[1], sum[0], 0, -c, 0, x);
            fmpz_neg_uiui(f, sum[1], sum[0]);
        }
        else                    /* coeff is non-negative, x non-negative */
        {
            if (x < c)
                fmpz_set_ui(f, c - x);  /* won't be negative and is smaller than c */
            else
                fmpz_neg_ui(f, x - c);  /* positive or zero */
        }
    }
    else
    {
        __mpz_struct * mc, * mf;
        mf = _fmpz_promote(f);    /* g is already large */
        mc = COEFF_TO_PTR(c);
        flint_mpz_sub_ui(mf, mc, x);
        _fmpz_demote_val(f);    /* cancellation may have occurred */
    }
}

void
sample_add_new(void * arg, ulong count)
{
    fmpz_t res, a;
    ulong ix, jx, b;
    int bits = *((int *) arg);

    FLINT_TEST_INIT(state);

    fmpz_init(res);
    fmpz_init(a);
   
    for (ix = 0; ix < 100 * count; ix++)
    {
        fmpz_randtest(a, state, bits);
        b = n_randtest(state);

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_add_ui(res, a, b);
        prof_stop();
    }

    flint_randclear(state);
    fmpz_clear(res);
    fmpz_clear(a);
}

void
sample_add_old(void * arg, ulong count)
{
    fmpz_t res, a;
    ulong ix, jx, b;
    int bits = *((int *) arg);

    FLINT_TEST_INIT(state);

    fmpz_init(res);
    fmpz_init(a);
   
    for (ix = 0; ix < 100 * count; ix++)
    {
        fmpz_randtest(a, state, bits);
        b = n_randtest(state);

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_add_ui_old(res, a, b);
        prof_stop();
    }

    flint_randclear(state);
    fmpz_clear(res);
    fmpz_clear(a);
}

void
sample_sub_new(void * arg, ulong count)
{
    fmpz_t res, a;
    ulong ix, jx, b;
    int bits = *((int *) arg);

    FLINT_TEST_INIT(state);

    fmpz_init(res);
    fmpz_init(a);
   
    for (ix = 0; ix < 100 * count; ix++)
    {
        fmpz_randtest(a, state, bits);
        b = n_randtest(state);

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_sub_ui(res, a, b);
        prof_stop();
    }

    flint_randclear(state);
    fmpz_clear(res);
    fmpz_clear(a);
}

void
sample_sub_old(void * arg, ulong count)
{
    fmpz_t res, a;
    ulong ix, jx, b;
    int bits = *((int *) arg);

    FLINT_TEST_INIT(state);

    fmpz_init(res);
    fmpz_init(a);
   
    for (ix = 0; ix < 100 * count; ix++)
    {
        fmpz_randtest(a, state, bits);
        b = n_randtest(state);

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_sub_ui_old(res, a, b);
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

    flint_printf("ADD\n");
    for (bits = 5; bits <= 150; bits += 5)
    {
        prof_repeat(&minnew, &maxnew, sample_add_new, &bits);
        prof_repeat(&minold, &maxold, sample_add_old, &bits);
        
        flint_printf("%d bits:      min %.2fx,    max %.2fx\n",
                bits, minold / minnew, maxold / maxnew);
    }

    flint_printf("\nSUB\n");
    for (bits = 5; bits <= 150; bits += 5)
    {
        prof_repeat(&minnew, &maxnew, sample_sub_new, &bits);
        prof_repeat(&minold, &maxold, sample_sub_old, &bits);
        
        flint_printf("%d bits:      min %.2fx,    max %.2fx\n",
                bits, minold / minnew, maxold / maxnew);
    }

    return 0;
}
