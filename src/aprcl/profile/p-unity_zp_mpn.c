/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "profiler.h"
#include "fmpz.h"
#include "gr.h"
#include "mpn_mod.h"
#include "aprcl.h"

int main(void)
{
    slong bits, pi;
    ulong ps[4] = {4, 9, 5, 11};
    flint_rand_t state;

    flint_rand_init(state);

    flint_printf("unity_zp_mpn_sqr, ns per call:\n");
    flint_printf("%12s", "bits \\ p^k");
    for (pi = 0; pi < 4; pi++)
        flint_printf(" %9wu sq %6wu ml", ps[pi], ps[pi]);
    flint_printf("\n");

    for (bits = 80; bits <= 320; bits += (bits < 128 ? 16 : 64))
    {
        fmpz_t n;
        gr_ctx_t ctx;

        fmpz_init(n);
        fmpz_randbits_unsigned(n, state, bits);
        fmpz_setbit(n, bits - 1);
        fmpz_setbit(n, 0);

        GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(ctx, n));

        flint_printf("%12wd", bits);

        for (pi = 0; pi < 4; pi++)
        {
            unity_zp_mpn f, g;
            slong i, d, reps;
            timeit_t timer;

            {
                ulong p = ps[pi], e = 1;
                if (p == 4) { p = 2; e = 2; }
                if (p == 9) { p = 3; e = 2; }
                unity_zp_mpn_init(f, p, e, ctx);
                unity_zp_mpn_init(g, p, e, ctx);
            }
            d = f->d;

            for (i = 0; i < d; i++)
            {
                fmpz_t c;
                fmpz_init(c);
                fmpz_randtest_mod(c, state, n);
                GR_MUST_SUCCEED(mpn_mod_set_fmpz(g->coeffs
                    + i * MPN_MOD_CTX_NLIMBS(ctx), c, ctx));
                fmpz_clear(c);
            }

            reps = 30000000 / (d * d * bits / 64) + 1;

            timeit_start(timer);
            for (i = 0; i < reps; i++)
                unity_zp_mpn_sqr(f, g);
            timeit_stop(timer);

            flint_printf(" %9.0f", timer->wall * 1e6 / reps);

            timeit_start(timer);
            for (i = 0; i < reps; i++)
                unity_zp_mpn_mul(f, g, g);
            timeit_stop(timer);

            flint_printf(" %9.0f", timer->wall * 1e6 / reps);

            unity_zp_mpn_clear(f);
            unity_zp_mpn_clear(g);
        }

        flint_printf("\n");
        gr_ctx_clear(ctx);
        fmpz_clear(n);
    }

    flint_rand_clear(state);
    return 0;
}
