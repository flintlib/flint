/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint/flint.h"
#include "flint/ulong_extras.h"
#include "flint/fmpz.h"
#include "flint/profiler.h"

typedef struct
{
   slong bits;
}
info_t;

ulong
z_gcd_old(slong a, slong b)
{
    ulong ua = FLINT_ABS(a);
    ulong ub = FLINT_ABS(b);

    return n_gcd(ua, ub);
}

void
fmpz_gcd_old(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (fmpz_is_zero(g))
    {
        fmpz_abs(f, h);
        return;
    }

    if (fmpz_is_zero(h))
    {
        fmpz_abs(f, g);
        return;
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (!COEFF_IS_MPZ(c2))  /* h is also small */
        {
            fmpz_set_si(f, z_gcd_old(c1, c2));
        }
        else                    /* h is large, but g is small */
        {
            fmpz c2d = fmpz_fdiv_ui(h, FLINT_ABS(c1));
            fmpz_set_si(f, z_gcd_old(c1, c2d));
        }
    }
    else
    {
        if (!COEFF_IS_MPZ(c2))  /* h is small, but g is large */
        {
            fmpz c1d = fmpz_fdiv_ui(g, FLINT_ABS(c2));
            fmpz_set_si(f, z_gcd_old(c2, c1d));
        }
        else                    /* g and h are both large */
        {
            __mpz_struct *mpz_ptr = _fmpz_promote(f);   /* aliasing fine as g, h already large */

            mpz_gcd(mpz_ptr, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
            _fmpz_demote_val(f);    /* gcd may be small */
        }
    }
}


void
sample_new(void * arg, ulong count)
{
    fmpz_t r, a, b;
    int ix;
    info_t * info = (info_t *) arg;
    slong bits = info->bits;

    fmpz_init(r);
    fmpz_init(a);
    fmpz_init(b);

    FLINT_TEST_INIT(state);

    prof_start();
    for (ix = 0; ix < 1000 * count; ix++)
    {
        fmpz_randtest(a, state, bits);
        fmpz_randtest(b, state, bits);
        fmpz_gcd(r, a, b);
    }
    prof_stop();

    fmpz_clear(r);
    fmpz_clear(a);
    fmpz_clear(b);

    flint_randclear(state);
}

void
sample_old(void * arg, ulong count)
{
    fmpz_t r, a, b;
    int ix;
    info_t * info = (info_t *) arg;
    slong bits = info->bits;

    fmpz_init(r);
    fmpz_init(a);
    fmpz_init(b);

    FLINT_TEST_INIT(state);

    prof_start();
    for (ix = 0; ix < 1000 * count; ix++)
    {
        fmpz_randtest(a, state, bits);
        fmpz_randtest(b, state, bits);
        fmpz_gcd_old(r, a, b);
    }
    prof_stop();

    fmpz_clear(r);
    fmpz_clear(a);
    fmpz_clear(b);

    flint_randclear(state);
}

int
main()
{
    double minnew, maxnew, minold, maxold;
    int bits;
    info_t as;

    for (bits = 5; bits <= 150; bits += 5)
    {
        as.bits = bits;

        prof_repeat(&minnew, &maxnew, sample_new, &as);
        prof_repeat(&minold, &maxold, sample_old, &as);
        flint_printf("%3d bits:   (min) %.2fx speedup   (max) %.2fx speedup\n",
                bits, minold / minnew, maxold / maxnew);
    }
}
