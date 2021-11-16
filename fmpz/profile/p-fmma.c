/*
    Copyright 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/profiler.h"

void
fmpz_fmma_old(fmpz_t f, const fmpz_t a, const fmpz_t b,
                        const fmpz_t c, const fmpz_t d)
{
    fmpz s, t, u, v;

    s = *a;
    t = *b;
    u = *c;
    v = *d;

    if (s == 0 || t == 0)
    {
        fmpz_mul(f, c, d);
        return;
    }

    if (u == 0 || v == 0)
    {
        fmpz_mul(f, a, b);
        return;
    }

    if (!COEFF_IS_MPZ(s) && !COEFF_IS_MPZ(t) &&
        !COEFF_IS_MPZ(u) && !COEFF_IS_MPZ(v))
    {
        mp_limb_t sh, sl, th, tl;

        smul_ppmm(sh, sl, s, t);
        smul_ppmm(th, tl, u, v);
        add_ssaaaa(sh, sl, sh, sl, th, tl);

        fmpz_set_signed_uiui(f, sh, sl);
        return;
    }

    if (f == c || f == d)
    {
        if (f == a || f == b)
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_mul(t, a, b);
            fmpz_addmul(t, c, d);
            fmpz_swap(t, f);
            fmpz_clear(t);
        }
        else
        {
            fmpz_mul(f, c, d);
            fmpz_addmul(f, a, b);
        }
    }
    else
    {
        fmpz_mul(f, a, b);
        fmpz_addmul(f, c, d);
    }
}

void sample_small(void * arg, ulong count)
{
    fmpz_t r, a, b, c, d;

    FLINT_TEST_INIT(state);

    fmpz_init(r);
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);
    fmpz_init(d);
    
    prof_start();
    for (int ix = 0; ix < 1000000 * count; ix++)
    {
        fmpz_set_ui(a, n_randint(state, COEFF_MAX));
        fmpz_set_ui(b, n_randint(state, COEFF_MAX));
        fmpz_set_ui(c, n_randint(state, COEFF_MAX));
        fmpz_set_ui(d, n_randint(state, COEFF_MAX));
        fmpz_fmma(r, a, b, c, d);
    }
    prof_stop();

    fmpz_clear(r);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
    fmpz_clear(d);

    flint_randclear(state);
}

void sample_small_old(void * arg, ulong count)
{
    fmpz_t r, a, b, c, d;

    FLINT_TEST_INIT(state);

    fmpz_init(r);
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);
    fmpz_init(d);
    
    prof_start();
    for (int ix = 0; ix < 1000000 * count; ix++)
    {
        fmpz_set_ui(a, n_randint(state, COEFF_MAX));
        fmpz_set_ui(b, n_randint(state, COEFF_MAX));
        fmpz_set_ui(c, n_randint(state, COEFF_MAX));
        fmpz_set_ui(d, n_randint(state, COEFF_MAX));
        fmpz_fmma_old(r, a, b, c, d);
    }
    prof_stop();

    fmpz_clear(r);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
    fmpz_clear(d);

    flint_randclear(state);
}

void sample_small_zeros(void * arg, ulong count)
{
    fmpz_t r, a, b, c, d;

    FLINT_TEST_INIT(state);

    fmpz_init(r);
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);
    fmpz_init(d);

    fmpz_zero(c);
    
    prof_start();
    for (int ix = 0; ix < 1000000 * count; ix++)
    {
        fmpz_set_ui(a, n_randint(state, COEFF_MAX));
        fmpz_set_ui(b, n_randint(state, COEFF_MAX));
        fmpz_set_ui(d, n_randint(state, COEFF_MAX));
        fmpz_fmma(r, a, b, c, d);
    }
    prof_stop();

    fmpz_clear(r);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
    fmpz_clear(d);

    flint_randclear(state);
}

void sample_small_zeros_old(void * arg, ulong count)
{
    fmpz_t r, a, b, c, d;

    FLINT_TEST_INIT(state);

    fmpz_init(r);
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);
    fmpz_init(d);

    fmpz_zero(c);
    
    prof_start();
    for (int ix = 0; ix < 1000000 * count; ix++)
    {
        fmpz_set_ui(a, n_randint(state, COEFF_MAX));
        fmpz_set_ui(b, n_randint(state, COEFF_MAX));
        fmpz_set_ui(d, n_randint(state, COEFF_MAX));
        fmpz_fmma_old(r, a, b, c, d);
    }
    prof_stop();

    fmpz_clear(r);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
    fmpz_clear(d);

    flint_randclear(state);
}

void sample_big_zeros(void * arg, ulong count)
{
    fmpz_t r, a, b, c, d;

    FLINT_TEST_INIT(state);

    fmpz_init(r);
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);
    fmpz_init(d);

    fmpz_zero(c);
    
    prof_start();
    for (int ix = 0; ix < 1000000 * count; ix++)
    {
        fmpz_randtest(a, state, 100);
        fmpz_randtest(b, state, 100);
        fmpz_randtest(d, state, 100);
        fmpz_fmma(r, a, b, c, d);
    }
    prof_stop();

    fmpz_clear(r);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
    fmpz_clear(d);

    flint_randclear(state);
}

void sample_big_zeros_old(void * arg, ulong count)
{
    fmpz_t r, a, b, c, d;

    FLINT_TEST_INIT(state);

    fmpz_init(r);
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);
    fmpz_init(d);

    fmpz_zero(c);
    
    prof_start();
    for (int ix = 0; ix < 1000000 * count; ix++)
    {
        fmpz_randtest(a, state, 100);
        fmpz_randtest(b, state, 100);
        fmpz_randtest(d, state, 100);
        fmpz_fmma_old(r, a, b, c, d);
    }
    prof_stop();

    fmpz_clear(r);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
    fmpz_clear(d);

    flint_randclear(state);
}

int main(void)
{
    double minnew, maxnew, minold, maxold;

    prof_repeat(&minnew, &maxnew, sample_small, NULL);
    prof_repeat(&minold, &maxold, sample_small_old, NULL);
    flint_printf("fmpz_fmma with small numbers:\n"
                 "      min %.3fx speedup,    max %.3fx speedup\n",
                 minold / minnew, maxold / maxnew);

    prof_repeat(&minnew, &maxnew, sample_small_zeros, NULL);
    prof_repeat(&minold, &maxold, sample_small_zeros_old, NULL);
    flint_printf("fmpz_fmma with one zero and small entries:\n"
                 "      min %.3fx speedup,    max %.3fx speedup\n",
                 minold / minnew, maxold / maxnew);

    prof_repeat(&minnew, &maxnew, sample_big_zeros, NULL);
    prof_repeat(&minold, &maxold, sample_big_zeros_old, NULL);
    flint_printf("fmpz_fmma with one zero and big entries:\n"
                 "      min %.3fx speedup,    max %.3fx speedup\n",
                 minold / minnew, maxold / maxnew);

    return 0;
}
