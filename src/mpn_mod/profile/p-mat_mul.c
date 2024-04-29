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
            if (__timer->cpu >= 30) \
                break; \
            __reps *= 10; \
        } \
    } while (0);
#endif

slong parameters[] = { 4, 5, 6, 8, 10, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 72, 80, 88, 96, 104, 112, 128, 144, 160, 192, 224, 256, 0 };

void
randvec(gr_ptr vec, flint_rand_t state, slong len, gr_ctx_t ctx)
{
    slong i;
    fmpz_t t;

    fmpz_init(t);
    for (i = 0; i < len; i++)
    {
        fmpz_randbits(t, state, MPN_MOD_CTX_MODULUS_BITS(ctx) + 10);
        GR_IGNORE(gr_set_fmpz(GR_ENTRY(vec, i, ctx->sizeof_elem), t, ctx));
    }

    fmpz_clear(t);
}

void
randmat(gr_mat_t mat, flint_rand_t state, gr_ctx_t ctx)
{
    slong i;

    for (i = 0; i < mat->r; i++)
        randvec(mat->rows[i], state, mat->c, ctx);
}

slong cutoff_bits[1000];
slong waksman_cutoffs[1000];
slong multi_mod_cutoffs[1000];
slong num_cutoffs = 0;

void
cutoffs(flint_rand_t state, gr_ctx_t ctx)
{
    gr_mat_t A, B, C;
    double t_classical, t_waksman, t_multi_mod, best, __;
    slong n, ni;
    char * beststr = "none";
    slong bits = MPN_MOD_CTX_MODULUS_BITS(ctx);

    int found_nonclassical = 0;

    cutoff_bits[num_cutoffs] = bits;

    for (n = 4; ; n = n < 100 ? n + 1 : (n < 600 ? n + 10 : n + 100))
    {
        gr_mat_init(A, n, n, ctx);
        gr_mat_init(B, n, n, ctx);
        gr_mat_init(C, n, n, ctx);

        randmat(A, state, ctx);
        randmat(B, state, ctx);

        TIMEIT_START
        GR_MUST_SUCCEED(gr_mat_mul_classical(C, A, B, ctx));
        TIMEIT_STOP_VALUES(__, t_classical)

        TIMEIT_START
        GR_MUST_SUCCEED(mpn_mod_mat_mul_waksman(C, A, B, ctx));
        TIMEIT_STOP_VALUES(__, t_waksman)

        TIMEIT_START
        GR_MUST_SUCCEED(mpn_mod_mat_mul_multi_mod(C, A, B, ctx));
        TIMEIT_STOP_VALUES(__, t_multi_mod)

        best = FLINT_MIN(t_classical, t_waksman);
        best = FLINT_MIN(best, t_multi_mod);

        if (best == t_multi_mod) beststr = "multi_mod";
        if (best == t_waksman) beststr = "waksman";
        if (best == t_classical) beststr = "classical";

        flint_printf("%wd   %wd    %g  %g  %g     %s\n", bits, n, t_classical, t_waksman, t_multi_mod, beststr);

        (void) __;

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(C, ctx);

        if (!found_nonclassical)
        {
            if (best != t_classical)
            {
                waksman_cutoffs[num_cutoffs] = n;
                found_nonclassical = 1;
            }
        }

        if (best != t_classical && best != t_waksman && t_multi_mod < t_classical && t_multi_mod < t_waksman)
        {
            multi_mod_cutoffs[num_cutoffs] = n;

            num_cutoffs++;
            break;
        }
    }

    for (ni = 0; ni < num_cutoffs; ni++)
        flint_printf("  {%wd, %wd},   /* bits = %wd */\n", waksman_cutoffs[ni], multi_mod_cutoffs[ni], cutoff_bits[ni]);
}

int main()
{
    fmpz_t p;
    gr_ctx_t ctx;
    flint_rand_t state;
    slong bits;

    flint_randinit(state);

    /* flint_set_num_threads(8); */

    for (bits = 80; bits <= 1024; bits += 16)
    {
        fmpz_init(p);
        fmpz_randprime(p, state, bits, 0);
        GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(ctx, p));
        cutoffs(state, ctx);
        gr_ctx_clear(ctx);
    }

    fmpz_clear(p);
    flint_randclear(state);
}
