/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "acb.h"
#include "acb_dirichlet.h"
#include "gr.h"
#include "dfloat.h"

/* The sieved power sum: the dfloat versions (plain midpoints with a
   priori bounds, and all balls) overlap acb_dirichlet_powsum_sieved at
   high precision, for s = sigma + i t with sigma in [-2, 3] (exact or
   a ball), t up to 10^16 (exact or a ball) and N up to a few 10^4, at
   every pair of precisions; the radius is tight for exact s.  The
   generic implementation is also exercised with arb balls for the
   terms and phases, with arf floats for the composites. */

TEST_FUNCTION_START(powsum, state)
{
    slong iter;
    static const int pairs[5][2] = {{1, 2}, {1, 3}, {2, 3}, {2, 4}, {3, 4}};

    DFLOAT_TEST_SKIP_IF_UNSUPPORTED(state);

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        acb_t s, r, z;
        ulong N, M;
        int k, ball, ok;
        double sig;

        acb_init(s);
        acb_init(r);
        acb_init(z);

        switch (n_randint(state, 4))
        {
            case 0: N = n_randint(state, 20); break;
            case 1: N = n_randint(state, 3000); break;
            case 2: N = 1 + n_randint(state, 30000); break;
            default: N = (UWORD(1) << n_randint(state, 15)) - n_randint(state, 2); break;
        }
        /* the table bound: all of it, or the two-phase version */
        M = n_randint(state, 2) ? 0 : n_randint(state, N + 1);
        sig = n_randint(state, 2) ? 0.5 : ((double) n_randint(state, 21) - 8) / 4.0;
        arb_set_d(acb_realref(s), sig);
        arb_set_ui(acb_imagref(s), n_randint(state, 1000));
        arb_mul_2exp_si(acb_imagref(s), acb_imagref(s), n_randint(state, 44));
        {
            arb_t u;
            arb_init(u);
            arb_set_ui(u, n_randtest(state));
            arb_mul_2exp_si(u, u, -64);
            arb_add(acb_imagref(s), acb_imagref(s), u, 200);
            if (n_randint(state, 2))
                arb_neg(acb_imagref(s), acb_imagref(s));
            arb_clear(u);
        }
        if (n_randint(state, 4) == 0)
            mag_set_ui_2exp_si(arb_radref(acb_imagref(s)), 1, -60 - (slong) n_randint(state, 100));
        if (n_randint(state, 4) == 0)
            mag_set_ui_2exp_si(arb_radref(acb_realref(s)), 1, -60 - (slong) n_randint(state, 100));

        acb_dirichlet_powsum_sieved(r, s, N, 1, 300);

        for (k = 0; k < 5; k++)
        {
            for (ball = 0; ball <= 1; ball++)
            {
                if (ball)
                    ok = _dfloat_powsum_sieved_ball(z, s, N, pairs[k][0], pairs[k][1], M);
                else
                    ok = _dfloat_powsum_sieved(z, s, N, pairs[k][0], pairs[k][1], M);

                if (!ok)
                    continue;

                if (!acb_overlaps(z, r))
                {
                    flint_printf("FAIL: overlap (T = %d, P = %d, ball = %d)\n", pairs[k][0], pairs[k][1], ball);
                    flint_printf("N = %wu, s = ", N); acb_printd(s, 30); flint_printf("\n");
                    flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
                    flint_printf("r = "); acb_printd(r, 30); flint_printf("\n");
                    flint_abort();
                }

                /* tightness for exact s with a moderate phase: about
                   53 T bits relative to the sum of |n^-s| */
                if (acb_is_exact(s) && N >= 2 && sig == 0.5 && arf_cmpabs_2exp_si(arb_midref(acb_imagref(s)), 30) < 0)
                {
                    mag_t m;
                    double w = 2.0 * sqrt((double) N) + 1.0;
                    mag_init(m);
                    mag_max(m, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));
                    if (mag_cmp_2exp_si(m, (slong) (-53 * pairs[k][0] + 24 + 2 * log2(w))) > 0)
                    {
                        flint_printf("FAIL: radius (T = %d, P = %d, ball = %d)\n", pairs[k][0], pairs[k][1], ball);
                        flint_printf("N = %wu, s = ", N); acb_printd(s, 30); flint_printf("\n");
                        flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
                        flint_abort();
                    }
                    mag_clear(m);
                }
            }
        }

        /* the automatic choice */
        if (dfloat_powsum_sieved(z, s, N, 60) && !acb_overlaps(z, r))
        {
            flint_printf("FAIL: dfloat_powsum_sieved\n");
            flint_abort();
        }

        /* other rings: arb terms and phases, arf composites */
        if (N <= 3000 && n_randint(state, 4) == 0)
        {
            gr_ctx_t T, P, F;
            slong tprec = 64 + n_randint(state, 100);
            int status;
            gr_ctx_init_real_arb(T, tprec);
            gr_ctx_init_real_arb(P, tprec + 80);
            gr_ctx_init_real_float_arf(F, tprec);
            for (ball = 0; ball <= 1; ball++)
            {
                if (ball)
                    status = gr_powsum_sieved(z, s, N, T, P, NULL, 0.0, 0.0, M);
                else
                    status = gr_powsum_sieved(z, s, N, T, P, F, ldexp(1.0, -tprec) * (1.0 + 0x1p-40), ldexp(1.0, -tprec) * (1.0 + 0x1p-40), M);
                if (status != GR_SUCCESS)
                {
                    flint_printf("FAIL: status (arb, ball = %d)\n", ball);
                    flint_printf("N = %wu, s = ", N); acb_printd(s, 30); flint_printf("\n");
                    flint_abort();
                }
                if (!acb_overlaps(z, r))
                {
                    flint_printf("FAIL: overlap (arb, ball = %d)\n", ball);
                    flint_printf("N = %wu, s = ", N); acb_printd(s, 30); flint_printf("\n");
                    flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
                    flint_printf("r = "); acb_printd(r, 30); flint_printf("\n");
                    flint_abort();
                }
                if (acb_is_exact(s) && N >= 2 && sig == 0.5 && arf_cmpabs_2exp_si(arb_midref(acb_imagref(s)), 30) < 0)
                {
                    mag_t m;
                    double w = 2.0 * sqrt((double) N) + 1.0;
                    mag_init(m);
                    mag_max(m, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));
                    if (mag_cmp_2exp_si(m, (slong) (-tprec + 24 + 2 * log2(w))) > 0)
                    {
                        flint_printf("FAIL: radius (arb, ball = %d)\n", ball);
                        flint_printf("N = %wu, s = ", N); acb_printd(s, 30); flint_printf("\n");
                        flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
                        flint_abort();
                    }
                    mag_clear(m);
                }
            }
            gr_ctx_clear(T);
            gr_ctx_clear(P);
            gr_ctx_clear(F);
        }

        acb_clear(s);
        acb_clear(r);
        acb_clear(z);
    }

    /* the threaded second phase (several chunks): the same result on
       1 and 3 threads, overlapping the full-table version */
    for (iter = 0; iter < 2 * flint_test_multiplier(); iter++)
    {
        acb_t s, z1, z3, zf;
        ulong N = 600000 + n_randint(state, 600000);
        slong nt = flint_get_num_threads();
        int k = n_randint(state, 5), ball = (pairs[k][0] == 3) || n_randint(state, 2), ok1, ok3, okf;

        acb_init(s);
        acb_init(z1);
        acb_init(z3);
        acb_init(zf);
        arb_set_d(acb_realref(s), n_randint(state, 2) ? 0.5 : ((double) n_randint(state, 13) - 4) / 4.0);
        arb_set_ui(acb_imagref(s), n_randtest(state));
        arb_mul_2exp_si(acb_imagref(s), acb_imagref(s), -64 + (slong) n_randint(state, 46));

        flint_set_num_threads(1);
        ok1 = ball ? _dfloat_powsum_sieved_ball(z1, s, N, pairs[k][0], pairs[k][1], 0)
                   : _dfloat_powsum_sieved(z1, s, N, pairs[k][0], pairs[k][1], 0);
        flint_set_num_threads(3);
        ok3 = ball ? _dfloat_powsum_sieved_ball(z3, s, N, pairs[k][0], pairs[k][1], 0)
                   : _dfloat_powsum_sieved(z3, s, N, pairs[k][0], pairs[k][1], 0);
        flint_set_num_threads(nt);
        okf = ball ? _dfloat_powsum_sieved_ball(zf, s, N, pairs[k][0], pairs[k][1], N)
                   : _dfloat_powsum_sieved(zf, s, N, pairs[k][0], pairs[k][1], N);

        if (ok1 != ok3 || ok1 != okf || (ok1 && (!acb_equal(z1, z3) || !acb_overlaps(z1, zf))))
        {
            flint_printf("FAIL: threads (T = %d, P = %d, ball = %d)\n", pairs[k][0], pairs[k][1], ball);
            flint_printf("N = %wu, s = ", N); acb_printd(s, 30); flint_printf("\n");
            flint_printf("z1 = "); acb_printd(z1, 30); flint_printf("\n");
            flint_printf("z3 = "); acb_printd(z3, 30); flint_printf("\n");
            flint_printf("zf = "); acb_printd(zf, 30); flint_printf("\n");
            flint_abort();
        }

        acb_clear(s);
        acb_clear(z1);
        acb_clear(z3);
        acb_clear(zf);
    }

    TEST_FUNCTION_END(state);
}
