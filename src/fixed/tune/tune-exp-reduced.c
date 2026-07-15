/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Tune the algorithm choice of fixed_exp_reduced over a grid of
   (r, wn): time algorithms 1..4 for every pair, print the winner
   grid, and suggest the terms-coordinate thresholds
   (terms = 64 wn / r) that the automatic selection uses:
   FIXED_EXP_BURST_TERMS between sinh and one burst, and
   FIXED_EXP_FULLBURST_TERMS between one burst and the full
   bit-burst.  The sinh/direct threshold is the established
   EXP_USE_SINH crossover.

   Usage: tune-exp-reduced [wnmax] */

#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "mpn_extras.h"
#include "fixed.h"

static double
time_alg(nn_ptr y, nn_srcptr t, slong wn, flint_bitcnt_t r, int alg)
{
    double best = 1e300;
    int pass;
    for (pass = 0; pass < 3; pass++)
    {
        timeit_t tm;
        slong reps = 1, k;
        for (;;)
        {
            timeit_start(tm);
            for (k = 0; k < reps; k++)
                fixed_exp_reduced(y, t, wn, r, alg);
            timeit_stop(tm);
            if (tm->wall >= 20 || reps >= 1000000)
                break;
            reps *= 4;
        }
        best = FLINT_MIN(best, tm->wall * 1e-3 / reps);
    }
    return best;
}

int main(int argc, char ** argv)
{
    slong wnmax = (argc > 1) ? atol(argv[1]) : 8192;
    slong wns[64], nw = 0, i, j;
    flint_bitcnt_t rs[] = {16, 32, 64, 128, 192, 256, 384, 512, 768,
        1024, 1536, 2048, 3072, 4096};
    slong nr = sizeof(rs) / sizeof(rs[0]);
    flint_rand_t state;

    flint_rand_init(state);

    for (i = 16; i <= wnmax; i *= 2)
        wns[nw++] = i;

    flint_printf("winner (alg 1 direct / 2 sinh / 3 burst+sinh / "
        "4 full burst); terms = 64 wn / r\n\n");
    flint_printf("      wn:");
    for (j = 0; j < nw; j++)
        flint_printf("%7wd", wns[j]);
    flint_printf("\n");

    for (i = 0; i < nr; i++)
    {
        flint_printf("r=%6wd:", (slong) rs[i]);
        for (j = 0; j < nw; j++)
        {
            slong wn = wns[j];
            nn_ptr t, y;
            double ta[5];
            int alg, bestalg = 1;
            slong k;

            if (FLINT_BITS * (ulong) wn <= rs[i])
            {
                flint_printf("      -");
                continue;
            }

            t = flint_malloc((wn + 2) * sizeof(ulong));
            y = flint_malloc((wn + 2) * sizeof(ulong));
            flint_mpn_urandomb(t, state, FLINT_BITS * wn);
            for (k = 0; k < (slong) (rs[i] / FLINT_BITS) && k < wn; k++)
                t[wn - 1 - k] = 0;
            if ((rs[i] % FLINT_BITS)
                && (slong) (rs[i] / FLINT_BITS) < wn)
                t[wn - 1 - rs[i] / FLINT_BITS] >>= rs[i] % FLINT_BITS;

            for (alg = 1; alg <= 4; alg++)
            {
                ta[alg] = time_alg(y, t, wn, rs[i], alg);
                if (ta[alg] < ta[bestalg])
                    bestalg = alg;
            }
            flint_printf("      %d", bestalg);
            flint_free(t); flint_free(y);
        }
        flint_printf("\n");
    }

    flint_printf("\npaste thresholds after inspecting the grid:\n"
        "  -D FIXED_EXP_BURST_TERMS=...  (alg 2 -> 3 boundary)\n"
        "  -D FIXED_EXP_FULLBURST_TERMS=...  (alg 3 -> 4 boundary)\n");
    return 0;
}
