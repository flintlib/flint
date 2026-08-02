/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Profile a fixed-point elementary function against arb across
   precisions:

       p-fixed FUNC [nmax]

   with FUNC one of exp, log1p, atan, sin_cos, tan, exp_notab,
   sin_cos_notab.  Sizes run
   n = 1..12, then in geometric steps of about 4/3 (12, 16, 21, 28,
   ...), up to nmax (default 4096) or until arb has been at least as
   fast twice in a row -- whichever comes first.  For each size the
   program prints n, the precision in bits and decimal digits, the
   reduction parameter that r = 0 selects, the time per call of arb
   and of the fixed function, and the speedup ratio.

   Each function is called once at the target size before timing, so
   that the shared angle/logarithm tables are built outside the
   measurement (precomputation is worth optimizing but is ignored
   here).  The timing loop cycles over an array of random inputs
   rather than repeating one x: the bitwise reductions branch on the
   bits of the argument, and a single input would hide the branch
   misprediction penalties a mixed workload pays.

   The notab entries profile the table-free functions against
   arb_exp / arb_sin_cos with emphasis on LARGE sizes: they never
   stop early when arb is ahead (the interesting regime is the
   asymptotic one; pass nmax explicitly, e.g.
   p-fixed exp_notab 1000000), and the number of distinct inputs
   scales down with n to keep the input array within memory. */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "profiler.h"
#include "flint.h"
#include "mpn_extras.h"
#include "fixed.h"
#include "arb.h"

#define NUM_X 64

typedef struct
{
    const char * name;
    int outputs;                       /* 1 or 2 */
    void (*fixed1)(nn_ptr, nn_srcptr, slong, int);
    void (*fixed2)(nn_ptr, nn_ptr, nn_srcptr, slong, int);
    int (*arbf1)(arb_t, const arb_t, slong);
    void (*arbf2)(arb_t, arb_t, const arb_t, slong);
    int (*default_r)(slong);
    int no_early_stop;                 /* keep going while arb wins
                                          (the notab large-n runs) */
} which_t;

/* arb_* have void results; adapt through thin shims so the table
   entries share a signature */
static int a_exp(arb_t y, const arb_t x, slong p) { arb_exp(y, x, p); return 0; }
static int a_log1p(arb_t y, const arb_t x, slong p) { arb_log1p(y, x, p); return 0; }
static int a_atan(arb_t y, const arb_t x, slong p) { arb_atan(y, x, p); return 0; }
static int a_tan(arb_t y, const arb_t x, slong p) { arb_tan(y, x, p); return 0; }
static void a_sin_cos(arb_t s, arb_t c, const arb_t x, slong p) { arb_sin_cos(s, c, x, p); }

/* the notab functions choose r internally; adapt to the shared
   (n, r) signature and report r = 0 */
static void p_exp_notab(nn_ptr y, nn_srcptr x, slong n, int r)
{ (void) r; fixed_exp_notab(y, x, n); }
static void p_sin_cos_notab(nn_ptr s, nn_ptr c, nn_srcptr x, slong n,
    int r)
{ (void) r; fixed_sin_cos_notab(s, c, x, n); }
static int p_notab_r(slong n) { (void) n; return 0; }

static const which_t table[] = {
    { "exp", 1, fixed_exp_bitwise_rs, NULL, a_exp, NULL,
      fixed_exp_bitwise_rs_default_r },
    { "log1p", 1, fixed_log1p_bitwise_rs, NULL, a_log1p, NULL,
      fixed_log1p_bitwise_rs_default_r },
    { "atan", 1, fixed_atan_bitwise_rs, NULL, a_atan, NULL,
      fixed_atan_bitwise_rs_default_r },
    { "tan", 1, fixed_tan_bitwise_rs, NULL, a_tan, NULL,
      fixed_trig_bitwise_rs_default_r },
    { "sin_cos", 2, NULL, fixed_sin_cos_bitwise_rs, NULL, a_sin_cos,
      fixed_trig_bitwise_rs_default_r },
    { "exp_notab", 1, p_exp_notab, NULL, a_exp, NULL,
      p_notab_r, 1 },
    { "sin_cos_notab", 2, NULL, p_sin_cos_notab, NULL, a_sin_cos,
      p_notab_r, 1 },
};

/* time one pass over the input array, repeated until the wall clock
   is trustworthy; returns nanoseconds per call, best of three */
static double
time_fixed(const which_t * w, nn_ptr y1, nn_ptr y2,
    ulong (*xs)[NUM_X], slong n, slong num_x)
{
    double best = 1e300;
    int pass;

    for (pass = 0; pass < 3; pass++)
    {
        slong reps = 1;
        timeit_t tm;

        while (1)
        {
            slong it, k;

            timeit_start(tm);
            for (it = 0; it < reps; it++)
                for (k = 0; k < num_x; k++)
                {
                    if (w->outputs == 1)
                        w->fixed1(y1, &xs[0][0] + (size_t) k * n, n, 0);
                    else
                        w->fixed2(y1, y2, &xs[0][0] + (size_t) k * n, n, 0);
                }
            timeit_stop(tm);
            if (tm->wall >= 30)
                break;
            reps *= 4;
        }

        best = FLINT_MIN(best, 1e6 * tm->wall / (reps * num_x));
    }

    return best;
}

static double
time_arb(const which_t * w, arb_ptr xa, slong n, slong num_x)
{
    double best = 1e300;
    slong prec = FLINT_BITS * n;
    int pass;
    arb_t v, v2;

    arb_init(v);
    arb_init(v2);

    for (pass = 0; pass < 3; pass++)
    {
        slong reps = 1;
        timeit_t tm;

        while (1)
        {
            slong it, k;

            timeit_start(tm);
            for (it = 0; it < reps; it++)
                for (k = 0; k < num_x; k++)
                {
                    if (w->outputs == 1)
                        w->arbf1(v, xa + k, prec);
                    else
                        w->arbf2(v, v2, xa + k, prec);
                }
            timeit_stop(tm);
            if (tm->wall >= 30)
                break;
            reps *= 4;
        }

        best = FLINT_MIN(best, 1e6 * tm->wall / (reps * num_x));
    }

    arb_clear(v);
    arb_clear(v2);
    return best;
}

int
main(int argc, char * argv[])
{
    const which_t * w = NULL;
    slong nmax = (argc > 2) ? atol(argv[2]) : 4096;
    slong n, i, k, arb_won = 0;
    flint_rand_t state;

    for (i = 0; argc > 1 && i < (slong) (sizeof(table) / sizeof(table[0])); i++)
        if (strcmp(argv[1], table[i].name) == 0)
            w = table + i;

    if (w == NULL)
    {
        flint_printf("usage: p-fixed exp|log1p|atan|sin_cos|tan"
            "|exp_notab|sin_cos_notab [nmax]\n");
        return 1;
    }

    flint_rand_init(state);

    flint_printf("%6s %8s %8s %5s %12s %12s %8s\n",
        "n", "bits", "digits", "r", "arb (ns)", "fixed (ns)", "ratio");

    for (n = (FLINT_BITS == 32 && w->outputs == 2) ? 2 : 1;
         n <= nmax;
         n = (n < 12) ? n + 1 : FLINT_MAX(n + 1, n * 4 / 3))
    {
        ulong (*xs)[NUM_X];
        nn_ptr y1, y2;
        arb_ptr xa;
        fmpz_t man, e;
        double tf, ta;
        int r0 = w->default_r(n);
        slong num_x = FLINT_MAX(2,
            FLINT_MIN((slong) NUM_X, 4000000 / n));

        xs = flint_malloc((size_t) num_x * n * sizeof(ulong));
        y1 = flint_malloc((n + 2) * sizeof(ulong));
        y2 = flint_malloc((n + 2) * sizeof(ulong));
        xa = _arb_vec_init(num_x);
        fmpz_init(man);
        fmpz_init(e);
        fmpz_set_si(e, -FLINT_BITS * n);

        for (k = 0; k < num_x; k++)
        {
            for (i = 0; i < n; i++)
                (&xs[0][0])[(size_t) k * n + i] = n_randlimb(state);
            fmpz_set_ui_array(man, &xs[0][0] + (size_t) k * n, n);
            arb_set_fmpz_2exp(xa + k, man, e);
        }

        /* warm-up: one call builds the angle/log tables for this
           size so that precomputation stays out of the measurement */
        if (w->outputs == 1)
            w->fixed1(y1, &xs[0][0], n, 0);
        else
            w->fixed2(y1, y2, &xs[0][0], n, 0);

        tf = time_fixed(w, y1, y2, xs, n, num_x);
        ta = time_arb(w, xa, n, num_x);

        flint_printf("%6wd %8wd %8.0f %5d %12.1f %12.1f %8.3f\n",
            n, FLINT_BITS * n, FLINT_BITS * n * 0.30103, r0, ta, tf,
            ta / tf);
        fflush(stdout);

        _arb_vec_clear(xa, num_x);
        flint_free(xs);
        flint_free(y1);
        flint_free(y2);
        fmpz_clear(man);
        fmpz_clear(e);

        if (!w->no_early_stop)
        {
            if (ta <= tf)
            {
                if (++arb_won >= 2)
                    break;
            }
            else
                arb_won = 0;
        }
    }

    flint_rand_clear(state);
    return 0;
}
