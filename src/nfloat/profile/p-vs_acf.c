/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz.h"
#include "gr.h"
#include "gr_special.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "acf.h"
#include "nfloat.h"
#include "profiler.h"
#include "double_extras.h"

#if 1
#undef TIMEIT_END_REPEAT
#define TIMEIT_END_REPEAT(__timer, __reps) \
            } \
            timeit_stop(__timer); \
            if (__timer->cpu >= 100) \
                break; \
            __reps *= 10; \
        } \
    } while (0);
#endif

int main()
{
    gr_ptr vec1, vec2, vec3;
    gr_ptr x, I;
    gr_ctx_t ctx;
    int which;
    slong i, n;
    slong prec;
    double __, t, acf_tadd = 0.0, acf_tmul = 0.0, acf_tsqr = 0.0, acf_tmul_scalar = 0.0, acf_taddmul_scalar = 0.0, acf_tsum = 0.0, acf_tprod = 0.0, acf_tdot = 0.0;
    double nfloat_tadd = 0.0, nfloat_tmul = 0.0, nfloat_tsqr = 0.0, nfloat_tmul_scalar = 0.0, nfloat_taddmul_scalar = 0.0, nfloat_tsum = 0.0, nfloat_tprod = 0.0, nfloat_tdot = 0.0;

    flint_printf("                   _gr_vec_add          _gr_vec_mul       _gr_vec_sqr       _gr_vec_mul_scalar  _gr_vec_addmul_scalar  _gr_vec_sum          _gr_vec_product      _gr_vec_dot\n");

    for (prec = 64; prec <= 4096; prec = prec < 256 ? prec + 64 : prec * 2)
    {
        flint_printf("prec = %wd\n", prec);

        for (n = 10; n <= 100; n *= 10)
        {
            for (which = 0; which < 2; which++)
            {
                flint_rand_t state;
                flint_rand_init(state);

                if (which == 0)
                    gr_ctx_init_complex_float_acf(ctx, prec);
                else
                    nfloat_complex_ctx_init(ctx, prec, 0);

                x = gr_heap_init(ctx);
                I = gr_heap_init(ctx);
                vec1 = gr_heap_init_vec(n, ctx);
                vec2 = gr_heap_init_vec(n, ctx);
                vec3 = gr_heap_init_vec(n, ctx);

                GR_MUST_SUCCEED(gr_i(I, ctx));

                for (i = 0; i < n; i++)
                {
                    gr_ptr v1 = GR_ENTRY(vec1, i, ctx->sizeof_elem);
                    gr_ptr v2 = GR_ENTRY(vec2, i, ctx->sizeof_elem);

                    GR_MUST_SUCCEED(gr_set_si(v1, 1 + n_randint(state, 1000), ctx));
                    if (n_randint(state, 2))
                        GR_MUST_SUCCEED(gr_neg(v1, v1, ctx));
                    GR_MUST_SUCCEED(gr_mul(v1, v1, I, ctx));
                    GR_MUST_SUCCEED(gr_add_si(v1, v1, 1 + n_randint(state, 1000), ctx));
                    if (n_randint(state, 2))
                        GR_MUST_SUCCEED(gr_neg(v1, v1, ctx));
                    GR_MUST_SUCCEED(gr_div_ui(v1, v1, 1 + n_randint(state, 1000), ctx));

                    GR_MUST_SUCCEED(gr_set_si(v2, 1 + n_randint(state, 1000), ctx));
                    if (n_randint(state, 2))
                        GR_MUST_SUCCEED(gr_neg(v2, v2, ctx));
                    GR_MUST_SUCCEED(gr_mul(v2, v2, I, ctx));
                    GR_MUST_SUCCEED(gr_add_si(v2, v2, 1 + n_randint(state, 1000), ctx));
                    if (n_randint(state, 2))
                        GR_MUST_SUCCEED(gr_neg(v2, v2, ctx));
                    GR_MUST_SUCCEED(gr_div_ui(v2, v2, 1 + n_randint(state, 1000), ctx));

/*
                    if (n_randint(state, 10) == 0)
                        GR_MUST_SUCCEED(gr_zero(v2, ctx));
*/
                }

                TIMEIT_START
                GR_MUST_SUCCEED(_gr_vec_add(vec3, vec1, vec2, n, ctx));
                TIMEIT_STOP_VALUES(__, t)
                if (which == 0) acf_tadd = t; else nfloat_tadd = t;
                (void) __;

                TIMEIT_START
                GR_MUST_SUCCEED(_gr_vec_mul(vec3, vec1, vec2, n, ctx));
                TIMEIT_STOP_VALUES(__, t)
                if (which == 0) acf_tmul = t; else nfloat_tmul = t;
                (void) __;

                TIMEIT_START
                GR_MUST_SUCCEED(_gr_vec_mul(vec3, vec1, vec1, n, ctx));
                TIMEIT_STOP_VALUES(__, t)
                if (which == 0) acf_tsqr = t; else nfloat_tsqr = t;
                (void) __;

                GR_MUST_SUCCEED(gr_pi(x, ctx));
                GR_MUST_SUCCEED(gr_add(x, x, I, ctx));
                GR_MUST_SUCCEED(gr_inv(x, x, ctx));

                TIMEIT_START
                GR_MUST_SUCCEED(_gr_vec_mul_scalar(vec3, vec1, n, x, ctx));
                TIMEIT_STOP_VALUES(__, t)
                if (which == 0) acf_tmul_scalar = t; else nfloat_tmul_scalar = t;
                (void) __;

                TIMEIT_START
                GR_MUST_SUCCEED(_gr_vec_addmul_scalar(vec3, vec1, n, x, ctx));
                TIMEIT_STOP_VALUES(__, t)
                if (which == 0) acf_taddmul_scalar = t; else nfloat_taddmul_scalar = t;
                (void) __;

                TIMEIT_START
                GR_MUST_SUCCEED(_gr_vec_sum(x, vec1, n, ctx));
                TIMEIT_STOP_VALUES(__, t)
                if (which == 0) acf_tsum = t; else nfloat_tsum = t;
                (void) __;

                TIMEIT_START
                GR_MUST_SUCCEED(_gr_vec_product(x, vec1, n, ctx));
                TIMEIT_STOP_VALUES(__, t)
                if (which == 0) acf_tprod = t; else nfloat_tprod = t;

                TIMEIT_START
                GR_MUST_SUCCEED(_gr_vec_dot(x, NULL, 0, vec1, vec2, n, ctx));
                TIMEIT_STOP_VALUES(__, t)
                if (which == 0) acf_tdot = t; else nfloat_tdot = t;

                gr_heap_clear(x, ctx);
                gr_heap_clear(I, ctx);
                gr_heap_clear_vec(vec1, n, ctx);
                gr_heap_clear_vec(vec2, n, ctx);
                gr_heap_clear_vec(vec3, n, ctx);

                flint_rand_clear(state);
            }

            flint_printf("n = %4wd   ", n);
            flint_printf("     %.3e (%.3fx)   %.3e (%.3fx)   %.3e (%.3fx)   %.3e (%.3fx)   %.3e (%.3fx)   %.3e (%.3fx)   %.3e (%.3fx)   %.3e (%.3fx)\n",
                nfloat_tadd, acf_tadd / nfloat_tadd,
                nfloat_tmul, acf_tmul / nfloat_tmul,
                nfloat_tsqr, acf_tsqr / nfloat_tsqr,
                nfloat_tmul_scalar, acf_tmul_scalar / nfloat_tmul_scalar,
                nfloat_taddmul_scalar, acf_taddmul_scalar / nfloat_taddmul_scalar,
                nfloat_tsum, acf_tsum / nfloat_tsum,
                nfloat_tprod, acf_tprod / nfloat_tprod,
                nfloat_tdot, acf_tdot / nfloat_tdot);
        }
    }

    return 0;
}
