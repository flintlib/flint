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
#include "gr_vec.h"
#include "gr_mat.h"
#include "arf.h"
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
    gr_ptr vec, vec2;
    gr_ptr x;
    gr_ctx_t ctx;
    int which;
    slong i, n;
    slong prec;
    double __, t, arf_tsum = 0.0, arf_tprod = 0.0, arf_tdot = 0.0;
    double nfloat_tsum = 0.0, nfloat_tprod = 0.0, nfloat_tdot = 0.0;

    flint_printf("                   _gr_vec_sum          _gr_vec_product      _gr_vec_dot\n");

    for (prec = 64; prec <= 2048; prec *= 2)
    {
        flint_printf("prec = %wd\n", prec);

        for (n = 3; n <= 100; n = (n == 3) ? 10 : n * 10)
        {
            for (which = 0; which < 2; which++)
            {
                if (which == 0)
                    gr_ctx_init_real_float_arf(ctx, prec);
                else
                    nfloat_ctx_init(ctx, prec, 0);

                x = gr_heap_init(ctx);
                vec = gr_heap_init_vec(n, ctx);
                vec2 = gr_heap_init_vec(n, ctx);

                for (i = 0; i < n; i++)
                {
                    GR_MUST_SUCCEED(gr_one(GR_ENTRY(vec, i, ctx->sizeof_elem), ctx));
                    if (i % 3 == 0)
                        GR_MUST_SUCCEED(gr_neg(GR_ENTRY(vec, i, ctx->sizeof_elem), GR_ENTRY(vec, i, ctx->sizeof_elem), ctx));
                    GR_MUST_SUCCEED(gr_div_ui(GR_ENTRY(vec, i, ctx->sizeof_elem), GR_ENTRY(vec, i, ctx->sizeof_elem), 2 * i + 3, ctx));
                }

                for (i = 0; i < n; i++)
                {
                    GR_MUST_SUCCEED(gr_one(GR_ENTRY(vec2, i, ctx->sizeof_elem), ctx));
                    if (i % 5 == 0)
                        GR_MUST_SUCCEED(gr_neg(GR_ENTRY(vec2, i, ctx->sizeof_elem), GR_ENTRY(vec2, i, ctx->sizeof_elem), ctx));
                    GR_MUST_SUCCEED(gr_div_ui(GR_ENTRY(vec2, i, ctx->sizeof_elem), GR_ENTRY(vec2, i, ctx->sizeof_elem), 2 * i + 5, ctx));
                }

                TIMEIT_START
                GR_MUST_SUCCEED(_gr_vec_sum(x, vec, n, ctx));
                TIMEIT_STOP_VALUES(__, t)
                if (which == 0) arf_tsum = t; else nfloat_tsum = t;
                (void) __;

                TIMEIT_START
                GR_MUST_SUCCEED(_gr_vec_product(x, vec, n, ctx));
                TIMEIT_STOP_VALUES(__, t)
                if (which == 0) arf_tprod = t; else nfloat_tprod = t;

                TIMEIT_START
                GR_MUST_SUCCEED(_gr_vec_dot(x, NULL, 0, vec, vec2, n, ctx));
                TIMEIT_STOP_VALUES(__, t)
                if (which == 0) arf_tdot = t; else nfloat_tdot = t;

                gr_heap_clear(x, ctx);
                gr_heap_clear_vec(vec, n, ctx);
            }

            flint_printf("n = %4wd   ", n);
            flint_printf("     %.3e (%.3fx)   %.3e (%.3fx)   %.3e (%.3fx)\n",
                nfloat_tsum, arf_tsum / nfloat_tsum, nfloat_tprod, arf_tprod / nfloat_tprod, nfloat_tdot, arf_tdot / nfloat_tdot);
        }
    }

    return 0;
}
