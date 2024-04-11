/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr.h"
#include "gr_poly.h"
#include "mpn_mod.h"
#include "profiler.h"

#if 1
#undef TIMEIT_END_REPEAT
#define TIMEIT_END_REPEAT(__timer, __reps) \
            } \
            timeit_stop(__timer); \
            if (__timer->cpu >= 10) \
                break; \
            __reps *= 10; \
        } \
    } while (0);
#endif

ulong parameter(slong i)
{
    if (i == 0)
        return 1;
    if (i % 2 == 1)
        return UWORD(1) << ((i + 1) / 2);
    return UWORD(3) << ((i - 2) / 2);
}

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

slong
find_karatsuba_cutoff(slong n, flint_rand_t state, gr_ctx_t ctx)
{
    slong i, m, found = -1;
    gr_ptr A, B, C, D;
    double t2, __, tbest;
    int squaring = 1;

    A = gr_heap_init_vec(n, ctx);
    B = gr_heap_init_vec(n, ctx);
    C = gr_heap_init_vec(2 * n - 1, ctx);
    D = gr_heap_init_vec(2 * n - 1, ctx);

    randvec(A, state, n, ctx);
    randvec(B, state, n, ctx);

    TIMEIT_START
    GR_MUST_SUCCEED(_mpn_mod_poly_mullow_karatsuba(C, A, n, squaring ? A : B, n, 2 * n - 1, 2, ctx));
    TIMEIT_STOP_VALUES(tbest, __)

    for (i = 2; ; i++)
    {
        m = parameter(i);

        if (m > n)
            break;

        TIMEIT_START
        GR_MUST_SUCCEED(_mpn_mod_poly_mullow_karatsuba(C, A, n, squaring ? A : B, n, 2 * n - 1, m, ctx));
        TIMEIT_STOP_VALUES(t2, __)
        (void) __;

        if (t2 < tbest)
        {
            tbest = t2;
            found = m;
        }
    }

    gr_heap_clear_vec(A, n, ctx);
    gr_heap_clear_vec(B, n, ctx);
    gr_heap_clear_vec(C, n, ctx);
    gr_heap_clear_vec(D, n, ctx);

    return found;
}

int main()
{
    fmpz_t p;
    gr_ctx_t ctx;
    flint_rand_t state;
    slong bits, m1, m2, m3;

    flint_randinit(state);

    for (bits = 65; bits <= 1024; bits++)
    {
        fmpz_init(p);
        fmpz_randprime(p, state, bits, 0);
        GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(ctx, p));

        m1 = find_karatsuba_cutoff(30, state, ctx);
        m2 = find_karatsuba_cutoff(75, state, ctx);
        m3 = find_karatsuba_cutoff(123, state, ctx);

        flint_printf("%wd:  %wd  %wd  %wd\n", bits, m1, m2, m3);

        gr_ctx_clear(ctx);
    }

    fmpz_clear(p);
    flint_randclear(state);
}
