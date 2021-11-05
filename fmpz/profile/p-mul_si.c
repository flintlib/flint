/*
    Copyright (C) 2021 Albin Ahlbäck, William Hart, Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "profiler.h"
#include "flint.h"
#include "ulong_extras.h"
#include "long_extras.h"
#include "fmpz.h"

#define ntests 100000

FMPZ_INLINE void
fmpz_set_signed_uiui_new(fmpz_t r, ulong hi, ulong lo)
{
    if (((slong) hi) < 0)
    {
        hi = -hi - (lo != 0);
        lo = -lo;
        fmpz_neg_uiui(r, hi, lo);
    }
    else
    {
        fmpz_set_uiui(r, hi, lo);
    }
}

void
fmpz_mul_si_new(fmpz_t f, const fmpz_t g, slong x)
{
    fmpz c2 = *g;

    if (!COEFF_IS_MPZ(c2)) /* c2 is small */
    {
        mp_limb_t prod[2];
        mp_limb_t mc2 = c2;
        mp_limb_t mx = x;

        /* limb by limb multiply (assembly for most CPU's) */
        smul_ppmm(prod[1], prod[0], mc2, mx);
        fmpz_set_signed_uiui(f, prod[1], prod[0]);
    }
    else                        /* c2 is large */
    {
        if (x == 0)
            fmpz_zero(f);
        else
        {
            __mpz_struct *mpz_ptr = _fmpz_promote(f);   /* ok without val as if aliased both are large */
            flint_mpz_mul_si(mpz_ptr, COEFF_TO_PTR(c2), x);
        }
    }
}

void
fmpz_mul_si_old(fmpz_t f, const fmpz_t g, slong x)
{
    fmpz c2 = *g;

    if (x == 0)
    {
        fmpz_zero(f);
        return;
    }
    else if (!COEFF_IS_MPZ(c2)) /* c2 is small */
    {
        mp_limb_t prod[2];
        mp_limb_t uc2 = FLINT_ABS(c2);
        mp_limb_t ux = FLINT_ABS(x);

        /* unsigned limb by limb multiply (assembly for most CPU's) */
        umul_ppmm(prod[1], prod[0], uc2, ux);
        if ((c2 ^ x) < WORD(0))
            fmpz_neg_uiui(f, prod[1], prod[0]);
        else
            fmpz_set_uiui(f, prod[1], prod[0]);
    }
    else                        /* c2 is large */
    {
        __mpz_struct *mpz_ptr = _fmpz_promote(f);   /* ok without val as if aliased both are large */
        flint_mpz_mul_si(mpz_ptr, COEFF_TO_PTR(c2), x);
    }
}

void
sample(void * arg, ulong count)
{
    fmpz_t res, a;
    slong b;
    ulong ix, jx;

    FLINT_TEST_INIT(state);

    fmpz_init(res);
    fmpz_init(a);
   
    for (ix = 0; ix < count; ix++)
    {
        fmpz_randtest(a, state, FLINT_BITS - 1);
        b = z_randint(state, WORD_MAX);

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_mul_si_new(res, a, b);
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
    slong b;
    ulong ix, jx;

    FLINT_TEST_INIT(state);

    fmpz_init(res);
    fmpz_init(a);
   
    for (ix = 0; ix < count; ix++)
    {
        fmpz_randtest(a, state, FLINT_BITS - 1);
        b = z_randint(state, WORD_MAX);

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_mul_si_old(res, a, b);
        prof_stop();
    }

    flint_randclear(state);
    fmpz_clear(res);
    fmpz_clear(a);
}

int
main(void)
{
    double min, max;

    prof_repeat(&min, &max, sample, NULL);
    flint_printf("fmpz_mul_si_new min time is %.3f cycles, max time is %.3f cycles\n", 
            min/FLINT_CLOCK_SCALE_FACTOR/ntests, max/FLINT_CLOCK_SCALE_FACTOR/ntests);

    prof_repeat(&min, &max, sample_old, NULL);
    flint_printf("fmpz_mul_si_old min time is %.3f cycles, max time is %.3f cycles\n", 
            min/FLINT_CLOCK_SCALE_FACTOR/ntests, max/FLINT_CLOCK_SCALE_FACTOR/ntests);

    return 0;
}
