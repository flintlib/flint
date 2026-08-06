/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "profiler.h"
#include "fmpz_poly.h"

int main(void)
{
    gr_poly_t f, g, u, v, c, d, e;
    gr_ctx_t R, Rx, Rxy, R1;
    fmpz_poly_t F;
    double t1, t2, FLINT_SET_BUT_UNUSED(tt);

    flint_rand_t state;
    flint_rand_init(state);

    slong N, M, i;

    gr_ctx_init_nmod(R1, 17);

    int ring;

    flint_printf("   N     M       fmpz  fmpz_poly      fmpzi       nmod        arb    fq_nmod         nf     gr_poly(nmod)\n");

    /* Outer length */
    for (N = 8; N <= 256; N = (N < 4 ? N + 1 : N * 2))
    {
        /* Inner length */
        for (M = 1; M <= 256; M *= 2)
        {
            flint_printf("%4wd  %4wd  ", N, M);
            fflush(stdout);

            for (ring = 0; ring < 8; ring++)
            {
                switch (ring)
                {
                    case 0:
                        gr_ctx_init_fmpz(R);
                        break;
                    case 1:
                        gr_ctx_init_fmpz_poly(R);
                        break;
                    case 2:
                        gr_ctx_init_fmpzi(R);
                        break;
                    case 3:
                        gr_ctx_init_nmod(R, 31337);
                        break;
                    case 4:
                        gr_ctx_init_real_arb(R, 512);
                        break;
                    case 5:
                        gr_ctx_init_fq_nmod(R, 17, 10, "a");
                        break;
                    case 6:
                        fmpz_poly_init(F);
                        fmpz_poly_set_coeff_si(F, 0, 1);
                        fmpz_poly_set_coeff_si(F, 1, 1);
                        fmpz_poly_set_coeff_si(F, 2, 1);
                        gr_ctx_init_nf_fmpz_poly(R, F);
                        fmpz_poly_clear(F);
                        break;
                    case 7:
                        gr_ctx_init_gr_poly(R, R1);
                        break;
                }

                gr_ctx_init_gr_poly(Rx, R);
                gr_ctx_init_gr_poly(Rxy, Rx);
                GR_MUST_SUCCEED(gr_ctx_set_gen_name(Rx, "x"));
                GR_MUST_SUCCEED(gr_ctx_set_gen_name(Rxy, "y"));

                gr_init(c, Rx);
                gr_init(d, Rx);
                gr_init(e, Rx);

                gr_init(f, Rxy);
                gr_init(g, Rxy);
                gr_init(u, Rxy);
                gr_init(v, Rxy);

                for (i = 0; i < N; i++)
                {
                    /* randtest polys are no good, try something more realistic */
                    if (gr_ctx_has_real_prec(R) == T_TRUE)
                    {
                        GR_MUST_SUCCEED(gr_set_str(c, "pi + 3*x", Rx));
                        GR_MUST_SUCCEED(gr_add_ui(c, c, i, Rx));
                        GR_MUST_SUCCEED(gr_poly_log_series(c, c, M, R));
                        GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(f, i, c, Rx));

                        GR_MUST_SUCCEED(gr_set_str(c, "(1/pi) + 5*x", Rx));
                        GR_MUST_SUCCEED(gr_add_ui(c, c, i, Rx));
                        GR_MUST_SUCCEED(gr_poly_log_series(c, c, M, R));
                        GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(g, i, c, Rx));
                    }
                    else
                    {
                        GR_MUST_SUCCEED(gr_poly_randtest(c, state, M, R));
                        GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(f, i, c, Rx));
                        GR_MUST_SUCCEED(gr_poly_randtest(c, state, M, R));
                        GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(g, i, c, Rx));
                    }
                }

                TIMEIT_START;
                GR_MUST_SUCCEED(gr_poly_mullow_classical(u, f, g, WORD_MAX, Rx));
                TIMEIT_STOP_VALUES(tt, t1);
                TIMEIT_START;
                GR_MUST_SUCCEED(gr_poly_mullow(v, f, g, WORD_MAX, Rx));
                TIMEIT_STOP_VALUES(tt, t2);

                flint_printf("%8.3fx  ", t1 / t2);
                fflush(stdout);

                if (gr_poly_equal(u, v, Rx) == T_FALSE)
                    flint_abort();

                gr_init(c, Rx);
                gr_init(d, Rx);
                gr_init(e, Rx);

                gr_clear(f, Rxy);
                gr_clear(g, Rxy);
                gr_clear(u, Rxy);
                gr_clear(v, Rxy);

                gr_ctx_clear(Rxy);
                gr_ctx_clear(Rx);
                gr_ctx_clear(R);
            }

            flint_printf("\n");
        }
    }

    flint_rand_clear(state);
    return 0;
}

