/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* The exponential of short dyadic arguments, for tuning the choice
   between the reductions and the direct bit-burst in exp.c:

       p-exp_window LIMBS HI LEN

   each a comma-separated list (LEN 0: all the remaining bits), for
   x = c 2^-(HI + LEN) with c a random integer of LEN bits (so x lies in
   [2^-(HI+1), 2^-HI)).  The columns are seconds per call of
   mp_real_exp_bits on x, of the table-free kernel (bit-burst) and of
   the diophantine kernel on x as a fraction, then of mp_real_exp_bits
   on -x and of the bit-burst on x followed by a reciprocal. */

#include <stdlib.h>
#include <string.h>
#include "flint.h"
#include "fmpz.h"
#include "profiler.h"
#include "mp_real.h"

static flint_rand_t st;
static mp_real_t X, Y, T1;
static nn_ptr V, Yv;
static slong N, SIGN;

static void run(int method)
{
    ulong err;
    if (method == 0)
        mp_real_exp_bits(Y, X, FLINT_BITS * N);
    else if (method == 1 || method == 2)
    {
        if (method == 1)
            _mp_real_exp_notab(Yv, &err, V, N);
        else
            _mp_real_exp_diophantine(Yv, &err, V, N);
        if (SIGN < 0)
        {
            /* the reciprocal */
            _mp_real_set_mpn_2exp(T1, Yv, N + 1, -FLINT_BITS * N);
            mp_real_set_ui(Y, 1);
            mp_real_div(Y, Y, T1, N);
        }
    }
}

static double timeit(int method)
{
    slong reps = 1, i;
    timeit_t t;
    for (;;)
    {
        timeit_start(t);
        for (i = 0; i < reps; i++)
            run(method);
        timeit_stop(t);
        if (t->wall >= 200 || reps >= 64)
            return 1e-3 * t->wall / reps;
        reps *= 2;
    }
}

int main(int argc, char ** argv)
{
    slong ns[16], nn = 0, hs[16], nh = 0, ls[16], nl = 0, i, j, k;
    char * s;
    fmpz_t f;
    char * lists[3] = { "1000,4000,15577", "0,16,64,256,1024", "64,1024,16384,0" };
    for (i = 0; i < 3; i++)
        if (argc > i + 1)
            lists[i] = argv[i + 1];
    for (s = lists[0]; *s && nn < 16; ) { ns[nn++] = strtol(s, &s, 10); if (*s) s++; }
    for (s = lists[1]; *s && nh < 16; ) { hs[nh++] = strtol(s, &s, 10); if (*s) s++; }
    for (s = lists[2]; *s && nl < 16; ) { ls[nl++] = strtol(s, &s, 10); if (*s) s++; }
    flint_rand_init(st);
    fmpz_init(f);
    mp_real_init(X); mp_real_init(Y); mp_real_init(T1);
    flint_printf("%8s %6s %8s | %9s %9s %9s | %9s %9s\n", "limbs", "hi", "len", "exp+", "notab+", "dioph+", "exp-", "notab-+div");
    for (i = 0; i < nn; i++)
    {
        N = ns[i];
        V = flint_malloc(sizeof(ulong) * N);
        Yv = flint_malloc(sizeof(ulong) * (N + 1));
        for (j = 0; j < nh; j++)
            for (k = 0; k < nl; k++)
            {
                slong hi = hs[j], len = ls[k] ? ls[k] : FLINT_BITS * N - hi;
                double t[5];
                if (hi + len > FLINT_BITS * N)
                    continue;
                /* x in [2^-(hi+1), 2^-hi) with len significant bits */
                fmpz_randbits(f, st, len);
                fmpz_abs(f, f);
                fmpz_setbit(f, len - 1);
                fmpz_setbit(f, 0);
                mp_real_set_fmpz(X, f);
                mp_real_mul_2exp_si(X, X, -(hi + len));
                {
                    ulong er;
                    _mp_real_get_fixed(V, &er, X, N);
                }
                SIGN = 1;
                mp_real_exp_bits(Y, X, FLINT_BITS * N);   /* warm caches */
                t[0] = timeit(0);
                t[1] = timeit(1);
                t[2] = (N <= 65536) ? timeit(2) : -1;
                SIGN = -1;
                mp_real_neg(X, X);
                mp_real_exp_bits(Y, X, FLINT_BITS * N);
                t[3] = timeit(0);
                t[4] = timeit(1);
                flint_printf("%8wd %6wd %8wd | %9.4f %9.4f %9.4f | %9.4f %9.4f\n", N, hi, len, t[0], t[1], t[2], t[3], t[4]);
                fflush(stdout);
            }
        flint_free(V);
        flint_free(Yv);
    }
    return 0;
}
