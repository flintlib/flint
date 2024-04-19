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
#include "gr_mat.h"
#include "mpn_mod.h"
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

void
randvec(gr_ptr vec, flint_rand_t state, slong len, slong bits, gr_ctx_t ctx)
{
    slong i;
    fmpz_t t;

    fmpz_init(t);
    for (i = 0; i < len; i++)
    {
        fmpz_randbits(t, state, bits + 10);
        GR_IGNORE(gr_set_fmpz(GR_ENTRY(vec, i, ctx->sizeof_elem), t, ctx));
    }

    fmpz_clear(t);
}

void
randmat(gr_mat_t mat, flint_rand_t state, slong bits, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < mat->r; i++)
        randvec(mat->rows[i], state, mat->c, bits, ctx);
}

int main()
{
    fmpz_t p;
    gr_ctx_t ctx, ctx2;
    flint_rand_t state;
    slong n, bits;
    double t1, t2, __;
    double speedup_add, speedup_mul, speedup_solve;
    slong i;
    slong bits_tab[] = { 80, 128, 180, 256, 512, 1024, };

    flint_randinit(state);

    gr_mat_t A, B, C;
    gr_mat_t A2, B2, C2;

    for (i = 0; i < 6; i++)
    {
        bits = bits_tab[i];

        fmpz_init(p);
        fmpz_randprime(p, state, bits, 0);

        GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(ctx, p));
        gr_ctx_init_fmpz_mod(ctx2, p);

        GR_MUST_SUCCEED(gr_ctx_set_is_field(ctx, T_TRUE));
        GR_MUST_SUCCEED(gr_ctx_set_is_field(ctx2, T_TRUE));

        flint_printf("  bits    n   add    mul    solve\n");

        for (n = 1; n <= 2048; n *= 2)
        {
            gr_mat_init(A, n, n, ctx);
            gr_mat_init(B, n, n, ctx);
            gr_mat_init(C, n, n, ctx);

            randmat(A, state, bits, ctx);
            randmat(B, state, bits, ctx);

            gr_mat_init(A2, n, n, ctx2);
            gr_mat_init(B2, n, n, ctx2);
            gr_mat_init(C2, n, n, ctx2);

            randmat(A2, state, bits, ctx2);
            randmat(B2, state, bits, ctx2);

            TIMEIT_START
            GR_MUST_SUCCEED(gr_mat_add(C, A, B, ctx));
            TIMEIT_STOP_VALUES(t1, __)
            /* don't time allocating fmpzs */
            GR_MUST_SUCCEED(gr_mat_add(C2, A2, B2, ctx2));
            TIMEIT_START
            GR_MUST_SUCCEED(gr_mat_add(C2, A2, B2, ctx2));
            TIMEIT_STOP_VALUES(t2, __)
            speedup_add = t2 / t1;

            TIMEIT_START
            GR_MUST_SUCCEED(gr_mat_mul(C, A, B, ctx));
            TIMEIT_STOP_VALUES(t1, __)
            TIMEIT_START
            GR_MUST_SUCCEED(gr_mat_mul(C2, A2, B2, ctx2));
            TIMEIT_STOP_VALUES(t2, __)

            speedup_mul = t2 / t1;

            gr_mat_clear(A, ctx);
            gr_mat_clear(B, ctx);
            gr_mat_clear(C, ctx);

            gr_mat_clear(A2, ctx2);
            gr_mat_clear(B2, ctx2);
            gr_mat_clear(C2, ctx2);


            gr_mat_init(A, n, n, ctx);
            gr_mat_init(B, n, 1, ctx);
            gr_mat_init(C, n, 1, ctx);

            randmat(A, state, bits, ctx);
            randmat(B, state, bits, ctx);

            gr_mat_init(A2, n, n, ctx2);
            gr_mat_init(B2, n, 1, ctx2);
            gr_mat_init(C2, n, 1, ctx2);

            randmat(A2, state, bits, ctx2);
            randmat(B2, state, bits, ctx2);

            TIMEIT_START
            GR_MUST_SUCCEED(gr_mat_nonsingular_solve(C, A, B, ctx));
            TIMEIT_STOP_VALUES(t1, __)
            TIMEIT_START
            GR_MUST_SUCCEED(gr_mat_nonsingular_solve(C2, A2, B2, ctx2));
            TIMEIT_STOP_VALUES(t2, __)
            speedup_solve = t2 / t1;

            gr_mat_clear(A, ctx);
            gr_mat_clear(B, ctx);
            gr_mat_clear(C, ctx);

            gr_mat_clear(A2, ctx2);
            gr_mat_clear(B2, ctx2);
            gr_mat_clear(C2, ctx2);

            flint_printf("%5wd %5wd   %.3f  %.3f  %.3f\n",
                bits, n, speedup_add, speedup_mul, speedup_solve);

            (void) __;
        }

        gr_ctx_clear(ctx);
        gr_ctx_clear(ctx2);
    }

    fmpz_clear(p);
    flint_randclear(state);
    return 0;
}
