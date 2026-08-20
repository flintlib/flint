/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "mpn_extras.h"
#include "profiler.h"

/* Compares the powm stages against each other and against mpz_powm,
   for random bases and for base 3 (the PRP-shaped workload). Reports
   milliseconds per call and the per-exponent-bit cost in units of one
   flint_mpn_sqr at the modulus size. */

static double
time_one(void (*f)(nn_ptr, nn_srcptr, nn_srcptr, mp_size_t,
                   nn_srcptr, mp_size_t),
         nn_ptr r, nn_srcptr b, nn_srcptr e, mp_size_t en,
         nn_srcptr m, mp_size_t mn, slong reps)
{
    timeit_t t;
    slong i;
    timeit_start(t);
    for (i = 0; i < reps; i++)
        f(r, b, e, en, m, mn);
    timeit_stop(t);
    return (double) t->cpu / reps;
}

int main(void)
{
    flint_rand_t state;
    mp_size_t sizes[] = { 400, 800, 1600, 3200, 6400, 12800, 25600 };
    slong nsizes = 7, si;

    flint_rand_init(state);

    flint_printf("powm, 256-bit exponent, times in ms "
                 "(sqr = one flint_mpn_sqr at mn):\n");
    flint_printf("%10s %8s %9s %9s %9s %9s %9s %11s\n",
                 "bits", "sqr", "basecase", "redc", "(unused)", "mpz",
                 "redc(b=3)", "ms/bit");

    for (si = 0; si < nsizes; si++)
    {
        mp_size_t mn = sizes[si], en = 4;
        nn_ptr m, b, e, r, X;
        mpz_t mz, bz, ez, rz;
        timeit_t ts;
        slong i, reps = FLINT_MAX(1, 4000 / mn), sreps;
        double t_sqr, t_bc, t_redc, t_fft, t_mpz, t_fft3;
        flint_bitcnt_t ebits;

        m = flint_malloc(mn * sizeof(ulong));
        b = flint_malloc(mn * sizeof(ulong));
        e = flint_malloc(en * sizeof(ulong));
        r = flint_malloc(mn * sizeof(ulong));
        X = flint_malloc(2 * mn * sizeof(ulong));

        flint_mpn_rrandom(m, state, mn);
        m[0] |= 1;
        m[mn - 1] |= UWORD(1) << 62;
        flint_mpn_rrandom(b, state, mn);
        mpn_tdiv_qr(X, b, 0, b, mn, m, mn);
        flint_mpn_rrandom(e, state, en);
        e[en - 1] |= UWORD(1) << 62;
        ebits = (en - 1) * FLINT_BITS + FLINT_BIT_COUNT(e[en - 1]);

        mpz_init(mz); mpz_init(bz); mpz_init(ez); mpz_init(rz);
        mpz_import(mz, mn, -1, sizeof(ulong), 0, 0, m);
        mpz_import(bz, mn, -1, sizeof(ulong), 0, 0, b);
        mpz_import(ez, en, -1, sizeof(ulong), 0, 0, e);

        sreps = FLINT_MAX(10, 400000 / mn);
        timeit_start(ts);
        for (i = 0; i < sreps; i++)
            flint_mpn_sqr(X, b, mn);
        timeit_stop(ts);
        t_sqr = (double) ts->cpu / sreps;

        t_bc = time_one(_flint_mpn_powm_basecase, r, b, e, en, m, mn, reps);
        t_redc = time_one(_flint_mpn_powm_redc, r, b, e, en, m, mn, reps);
        t_fft = 0.0;
        {
            timeit_t t;
            timeit_start(t);
            for (i = 0; i < reps; i++)
                mpz_powm(rz, bz, ez, mz);
            timeit_stop(t);
            t_mpz = (double) t->cpu / reps;
        }

        flint_mpn_zero(b, mn);
        b[0] = 3;
        t_fft3 = time_one(_flint_mpn_powm_redc, r, b, e, en, m, mn, reps);

        flint_printf("%10wd %8.2f %9.1f %9.1f %9.1f %9.1f %9.1f %11.4f\n",
                     (slong) mn * FLINT_BITS, t_sqr, t_bc, t_redc, t_fft,
                     t_mpz, t_fft3, t_fft3 / ebits);

        mpz_clear(mz); mpz_clear(bz); mpz_clear(ez); mpz_clear(rz);
        flint_free(m); flint_free(b); flint_free(e);
        flint_free(r); flint_free(X);
    }

    flint_rand_clear(state);
    return 0;
}
