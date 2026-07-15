/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Choose the reduction parameter r as a function of the number of
   intended evaluations, at a fixed precision.

   The bitwise functions amortize an r-dependent precomputation: the
   angle/log tables hold r + 1 entries of n + 1 limbs, and building
   them costs the more the larger r is, while each evaluation costs
   the less (fewer series terms).  tune-bitwise-r picks the r that
   minimizes the pure evaluation time -- the right choice when the
   tables are reused many times.  For a caller who will evaluate the
   function N times at precision n, the right objective is instead

       total(r) = P(r) + N E(r),

   with P the precomputation and E the per-evaluation time.  This
   program measures P(r) and E(r) for every ladder parameter
   r in {64, 128, ..., rmax} at the given n (resetting the cached
   tables between candidates so that P is measured cold), discards
   dominated candidates, and prints the lower envelope of the lines
   total_r(N): the breakpoints N* = (P2 - P1) / (E1 - E2) at which
   the next larger r starts paying for its precomputation, both as a
   human-readable table and as a paste-able C initializer

       static const struct { double n_evals; int r; } ...[] =
           { { 1.0, 64 }, { 12.3, 192 }, ... };

   where each entry's r applies from n_evals evaluations upward.

   Usage: tune-r-amortized <exp|log1p|atan|sincos|tan> <n> [rmax]
   (rmax defaults to 768 and is clamped to FLINT_BITS n - 16). */

#include <stdlib.h>
#include <string.h>
#include "profiler.h"
#include "flint.h"
#include "mpn_extras.h"
#include "fixed.h"

void _fixed_exp_logs_ensure(slong nv, slong rc);
void _fixed_atans_ensure(slong nv, slong rc);

typedef enum { F_EXP, F_LOG1P, F_ATAN, F_SINCOS, F_TAN } which_t;

static void
tables_reset(which_t w)
{
    if (w == F_EXP || w == F_LOG1P)
    {
        flint_free(_fixed_exp_logs);
        _fixed_exp_logs = NULL;
        _fixed_exp_logs_n = 0;
        _fixed_exp_logs_r = 0;
    }
    else
    {
        flint_free(_fixed_atans);
        _fixed_atans = NULL;
        _fixed_atans_n = 0;
        _fixed_atans_r = 0;
    }
}

static void
eval_once(which_t w, nn_ptr y, nn_ptr y2, nn_srcptr x, slong n, int r)
{
    switch (w)
    {
        case F_EXP:    fixed_exp_bitwise_rs(y, x, n, r); break;
        case F_LOG1P:  fixed_log1p_bitwise_rs(y, x, n, r); break;
        case F_ATAN:   fixed_atan_bitwise_rs(y, x, n, r); break;
        case F_SINCOS: fixed_sin_cos_bitwise_rs(y, y2, x, n, r); break;
        case F_TAN:    fixed_tan_bitwise_rs(y, x, n, r); break;
    }
}

/* seconds for one precomputation at (n, r), tables reset first */
static double
time_precomp(which_t w, slong n, int r)
{
    double best = 1e300;
    int pass, passes = (n <= 2000) ? 2 : 1;
    for (pass = 0; pass < passes; pass++)
    {
        timeit_t tm;
        tables_reset(w);
        timeit_start(tm);
        if (w == F_EXP || w == F_LOG1P)
            _fixed_exp_logs_ensure(n, r);
        else
            _fixed_atans_ensure(n, r);
        timeit_stop(tm);
        best = FLINT_MIN(best, tm->wall * 1e-3);
    }
    return best;
}

/* seconds per evaluation, tables assumed hot for (n, r) */
static double
time_eval(which_t w, nn_ptr y, nn_ptr y2, nn_srcptr x, slong n, int r)
{
    double best = 1e300;
    int pass;

    eval_once(w, y, y2, x, n, r);   /* safety: tables surely hot */

    for (pass = 0; pass < 3; pass++)
    {
        timeit_t tm;
        slong reps = 1, k;
        double dt;

        for (;;)
        {
            timeit_start(tm);
            for (k = 0; k < reps; k++)
                eval_once(w, y, y2, x, n, r);
            timeit_stop(tm);
            if (tm->wall >= 50 || reps >= 1000000)
                break;
            reps *= 4;
        }
        dt = tm->wall * 1e-3 / reps;
        best = FLINT_MIN(best, dt);
    }
    return best;
}

int main(int argc, char ** argv)
{
    which_t w;
    const char * names[] = {"exp", "log1p", "atan", "sincos", "tan"};
    slong n, rmax = 768, nr = 0, i, j;
    int rs[64];
    double P[64], E[64];
    int keep[64];
    nn_ptr x, y, y2;
    flint_rand_t state;

    if (argc < 3)
    {
        flint_printf("usage: tune-r-amortized "
            "<exp|log1p|atan|sincos|tan> <n> [rmax]\n");
        return 1;
    }
    for (i = 0; i < 5; i++)
        if (!strcmp(argv[1], names[i]))
            break;
    if (i == 5)
    {
        flint_printf("unknown function %s\n", argv[1]);
        return 1;
    }
    w = (which_t) i;
    n = atol(argv[2]);
    if (argc > 3)
        rmax = atol(argv[3]);
    rmax = FLINT_MIN(rmax, FLINT_BITS * n - 16);

    flint_rand_init(state);
    x = flint_malloc(n * sizeof(ulong));
    y = flint_malloc((n + 2) * sizeof(ulong));
    y2 = flint_malloc((n + 2) * sizeof(ulong));
    flint_mpn_urandomb(x, state, FLINT_BITS * n);
    x[n - 1] |= (UWORD(1) << (FLINT_BITS - 2));
    x[n - 1] &= ~(UWORD(1) << (FLINT_BITS - 1));

    flint_printf("fixed_%s_bitwise_rs, n = %wd\n\n", names[w], n);
    flint_printf("    r    precomp(us)      eval(us)\n");

    for (i = 0; 64 * (i + 1) <= rmax; i++)
    {
        rs[nr] = 64 * (int) (i + 1);
        P[nr] = time_precomp(w, n, rs[nr]);
        E[nr] = time_eval(w, y, y2, x, n, rs[nr]);
        flint_printf("%5d %14.0f %13.1f\n",
            rs[nr], P[nr] * 1e6, E[nr] * 1e6);
        nr++;
    }

    /* drop dominated candidates (both P and E no better than some
       other candidate's) */
    for (i = 0; i < nr; i++)
    {
        keep[i] = 1;
        for (j = 0; j < nr; j++)
            if (j != i && P[j] <= P[i] && E[j] <= E[i]
                && (P[j] < P[i] || E[j] < E[i]))
                keep[i] = 0;
    }

    /* lower envelope of total_r(N) = P + N E over N >= 1, walking
       from the smallest-total-at-N=1 candidate toward smaller E */
    {
        double breaks[64];
        int    er[64];
        slong  ne = 0;
        double Ncur = 1.0;
        int    cur = -1;

        /* best candidate at N = 1 */
        {
            double bestv = 1e300;
            for (i = 0; i < nr; i++)
                if (keep[i] && P[i] + E[i] < bestv)
                {
                    bestv = P[i] + E[i];
                    cur = (int) i;
                }
        }
        breaks[ne] = 1.0;
        er[ne] = rs[cur];
        ne++;

        /* walk the envelope: from the current line, the earliest
           crossing against a strictly smaller-E candidate is the
           next envelope member */
        while (1)
        {
            double bestN = 1e300;
            int nexti = -1;

            for (i = 0; i < nr; i++)
                if (keep[i] && E[i] < E[cur])
                {
                    double Nc = (P[i] - P[cur]) / (E[cur] - E[i]);
                    if (Nc >= Ncur && Nc < bestN)
                    {
                        bestN = Nc;
                        nexti = (int) i;
                    }
                }
            if (nexti < 0)
                break;
            cur = nexti;
            Ncur = bestN;
            breaks[ne] = Ncur;
            er[ne] = rs[cur];
            ne++;
        }

        flint_printf("\noptimal r by evaluation count:\n");
        for (i = 0; i < ne; i++)
        {
            if (i + 1 < ne)
                flint_printf("    %.4g <= N < %.4g evals: r = %d\n",
                    breaks[i], breaks[i + 1], er[i]);
            else
                flint_printf("    N >= %.4g evals: r = %d\n",
                    breaks[i], er[i]);
        }

        flint_printf("\n/* fixed_%s_bitwise_rs, n = %wd: smallest "
            "total-time r for N\n   evaluations (P + N E, measured); "
            "each entry applies from\n   n_evals upward */\n",
            names[w], n);
        flint_printf("static const struct { double n_evals; int r; }\n"
            "    %s_r_amortized_%wd[] = {\n    ", names[w], n);
        for (i = 0; i < ne; i++)
            flint_printf("{ %.4g, %d }, ", breaks[i], er[i]);
        flint_printf("\n};\n");
    }

    flint_free(x);
    flint_free(y);
    flint_free(y2);
    flint_rand_clear(state);
    return 0;
}
