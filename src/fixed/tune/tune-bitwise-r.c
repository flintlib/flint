/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Tune the LARGE-n default reduction parameters of the bitwise
   functions (used when they are called with r = 0 beyond the sizes
   the fully specialized *_opt_<n>.c files cover; those small-n r are
   tuned per size, at arbitrary r, by dev/tune_fixed.py -- this
   program picks from the ladder where the general series evaluation
   is available in the library).

   For each function, the candidate reduction parameters form a
   ladder (32, then every multiple of 64 up to rmax), and for each
   consecutive pair the program finds the smallest n at which the
   larger parameter stops losing, by binary search under the
   assumption that the optimal r is nondecreasing in n (which holds
   to within noise: the margins near the crossovers are flat).  A
   candidate only counts once it is genuinely available, i.e. not
   clamped by r <= FLINT_BITS n - 16.

   fixed_sin_cos_bitwise_rs and fixed_tan_bitwise_rs share the
   half-angle reduction and reconstruction, hence a single table
   (_fixed_trig_bitwise_rs_r_tab); it is tuned on sin_cos, the
   heavier and more common call, and the tangent run below is
   printed as a cross-check only -- its crossovers should agree to
   within noise.

   Usage: tune-bitwise-r [nmax] [rmax]     (defaults 640, 768;
   nmax at most 32767 so that the crossovers fit the short tables)

   The output is a pair of tables per function in a form that can be
   pasted directly into src/fixed/exp_bitwise_rs.c,
   src/fixed/log1p_bitwise_rs.c, src/fixed/atan_bitwise_rs.c and
   src/fixed/tan_bitwise_rs.c. */

#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "mpn_extras.h"
#include "fixed.h"

/* nanoseconds per call, minimum of npasses repetition-scaled runs */
static double
time_call(void (*func)(nn_ptr, nn_srcptr, slong, int), nn_ptr y,
          nn_srcptr x, slong n, int r, int npasses)
{
    double best = 1e300;
    int pass;

    /* one untimed call builds the angle/log tables for this (n, r)
       so that precomputation never leaks into the measurement */
    func(y, x, n, r);

    for (pass = 0; pass < npasses; pass++)
    {
        slong reps = 1;
        timeit_t tm;

        while (1)
        {
            slong ix;

            timeit_start(tm);
            for (ix = 0; ix < reps; ix++)
                func(y, x, n, r);
            timeit_stop(tm);
            if (tm->wall >= 5)
                break;
            reps *= 4;
        }

        best = FLINT_MIN(best, 1e6 * tm->wall / reps);
    }

    return best;
}

/* fixed_sin_cos_bitwise_rs computes two outputs; adapt it to the
   common signature by writing the cosine to a scratch buffer (both
   outputs are always computed internally, so this measures the full
   cost).  On 32-bit limbs n = 1 is unsupported: the probe is
   skipped there, which only affects the leading table entry, whose
   value is the r floor by convention. */
static nn_ptr sincos_scratch = NULL;
static slong sincos_scratch_n = 0;

static void
sin_cos_wrapper(nn_ptr y, nn_srcptr x, slong n, int r)
{
    if (FLINT_BITS == 32 && n < 2)
        return;

    if (n + 1 > sincos_scratch_n)
    {
        flint_free(sincos_scratch);
        sincos_scratch = flint_malloc((n + 1) * sizeof(ulong));
        sincos_scratch_n = n + 1;
    }

    fixed_sin_cos_bitwise_rs(y, sincos_scratch, x, n, r);
}

/* does r_new run at least about as fast as r_old at size n?
   Deterministic at the degenerate end: when r_new exceeds
   FLINT_BITS n - 16 it is not genuinely available and loses. */
static int
wins(void (*func)(nn_ptr, nn_srcptr, slong, int), nn_ptr y,
     nn_srcptr x, slong n, int r_new, int r_old)
{
    double t_new, t_old;

    if ((slong) r_new > FLINT_BITS * n - 16)
        return 0;

    t_new = time_call(func, y, x, n, r_new, 3);
    t_old = time_call(func, y, x, n, r_old, 3);

    return t_new <= 1.01 * t_old;
}

static void
tune_func(const char * name,
          void (*func)(nn_ptr, nn_srcptr, slong, int),
          slong nmax, slong rmax)
{
    slong maxtab = rmax / FLINT_BITS + 3;
    int * r_tab;
    slong * n_tab;
    slong num, j, k;
    nn_ptr x, y;
    flint_rand_t state;

    r_tab = flint_malloc(maxtab * sizeof(int));
    n_tab = flint_malloc(maxtab * sizeof(slong));
    x = flint_malloc(nmax * sizeof(ulong));
    y = flint_malloc((nmax + 1) * sizeof(ulong));

    flint_rand_init(state);
    /* uniform random limbs: structured random (bit runs) would
       degenerate the greedy accept pattern and skew the reduction
       cost */
    for (j = 0; j < nmax; j++)
        x[j] = n_randlimb(state);

    r_tab[0] = 32;
    n_tab[0] = 1;
    num = 1;

    while (num < maxtab)
    {
        int r_old = r_tab[num - 1];
        int r_new = (r_old < 64) ? 64 : r_old + 64;
        slong lo, hi;

        if ((slong) r_new > rmax)
            break;

        /* smallest n in (n_tab[num-1], nmax] where r_new wins;
           the optimum is assumed nondecreasing in n */
        if (!wins(func, y, x, nmax, r_new, r_old))
            break;

        lo = n_tab[num - 1];        /* r_new loses (or ties clamped) */
        hi = nmax;                  /* r_new wins */
        while (hi - lo > 1)
        {
            slong mid = lo + (hi - lo) / 2;

            if (wins(func, y, x, mid, r_new, r_old))
                hi = mid;
            else
                lo = mid;
        }

        flint_printf("/* %s: r = %d from n = %wd */\n",
            name, r_new, hi);
        fflush(stdout);

        r_tab[num] = r_new;
        n_tab[num] = hi;
        num++;
    }

    flint_printf("\n/* smallest n at which each reduction parameter"
        " becomes optimal\n   (generated by"
        " src/fixed/tune/tune-bitwise-r.c) */\n");
    flint_printf("static const int _%s_r_tab[] = {", name);
    for (j = 0; j < num; j++)
        flint_printf("%s%d", j ? ", " : "", r_tab[j]);
    flint_printf("};\n");
    flint_printf("static const short _%s_n_tab[] = {", name);
    for (j = 0; j < num; j++)
        flint_printf("%s%wd", j ? ", " : "", n_tab[j]);
    flint_printf("};\n\n");

    /* sanity sweep: print the selected r over a few sizes */
    flint_printf("/* selection:");
    for (k = 1; k <= nmax; k = (k < 8) ? k + 1 : k * 3 / 2)
    {
        for (j = 0; j + 1 < num && k >= n_tab[j + 1]; j++)
            ;
        flint_printf(" n=%wd:r=%d", k, r_tab[j]);
    }
    flint_printf(" */\n\n");

    flint_rand_clear(state);
    flint_free(r_tab);
    flint_free(n_tab);
    flint_free(x);
    flint_free(y);
}

int
main(int argc, char * argv[])
{
    slong nmax = (argc > 1) ? atol(argv[1]) : 640;
    slong rmax = (argc > 2) ? atol(argv[2]) : 768;

    if (nmax < 2 || nmax > 32767 || rmax < 32)
    {
        flint_printf("usage: tune-bitwise-r [nmax (2..32767)]"
            " [rmax (>= 32)]\n");
        return 1;
    }

    flint_printf("/* tuning the bitwise default r"
        "\n   (nmax = %wd, rmax = %wd) */\n\n", nmax, rmax);

    tune_func("fixed_exp_bitwise_rs", fixed_exp_bitwise_rs,
        nmax, rmax);
    tune_func("fixed_log1p_bitwise_rs", fixed_log1p_bitwise_rs,
        nmax, rmax);
    tune_func("fixed_atan_bitwise_rs", fixed_atan_bitwise_rs,
        nmax, rmax);
    /* sin_cos and tan share the half-angle path: the sin_cos run
       below yields _fixed_trig_bitwise_rs_{r,n}_tab, and the tangent
       run is a cross-check whose crossovers should agree */
    tune_func("fixed_trig_bitwise_rs", sin_cos_wrapper,
        nmax, rmax);
    tune_func("fixed_tan_bitwise_rs_CROSSCHECK", fixed_tan_bitwise_rs,
        nmax, rmax);

    flint_free(sincos_scratch);
    sincos_scratch = NULL;
    sincos_scratch_n = 0;

    return 0;
}
