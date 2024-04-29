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
#include "gr_poly.h"
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

int main()
{
    fmpz_t p;
    gr_ctx_t ctx, ctx2;
    flint_rand_t state;
    slong n, bits;
    slong lenA, lenB, lenC;
    double t1, t2, __;
    double speedup_mul, speedup_sqr, speedup_mullow, speedup_sqrlow, speedup_divrem, speedup_gcd;
    slong i;
    slong bits_tab[] = { 80, 128, 180, 256, 512, 1024, };

    flint_randinit(state);

    gr_ptr A, B, C, Q, R;
    gr_ptr A2, B2, C2, Q2, R2;

    for (i = 0; i < 6; i++)
    {
        bits = bits_tab[i];

        fmpz_init(p);
        fmpz_randprime(p, state, bits, 0);

        GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(ctx, p));
        gr_ctx_init_fmpz_mod(ctx2, p);

        flint_printf("  bits    n   mul    sqr    mullow sqrlow divrem gcd\n");

        for (n = 1; n <= 8192; n *= 2)
        {
            lenA = lenB = n;
            lenC = lenA + lenB - 1;

            A = gr_heap_init_vec(lenA, ctx);
            B = gr_heap_init_vec(lenB, ctx);
            C = gr_heap_init_vec(lenC, ctx);
            Q = gr_heap_init_vec(lenA, ctx);
            R = gr_heap_init_vec(lenB, ctx);

            randvec(A, state, lenA, bits, ctx);
            randvec(B, state, lenB, bits, ctx);

            A2 = gr_heap_init_vec(lenA, ctx2);
            B2 = gr_heap_init_vec(lenB, ctx2);
            C2 = gr_heap_init_vec(lenC, ctx2);
            Q2 = gr_heap_init_vec(lenA, ctx2);
            R2 = gr_heap_init_vec(lenB, ctx2);

            randvec(A2, state, lenA, bits, ctx2);
            randvec(B2, state, lenB, bits, ctx2);

            GR_MUST_SUCCEED(_gr_poly_mullow(C, A, lenA, B, lenB, lenA + lenB - 1, ctx));
            GR_MUST_SUCCEED(_gr_poly_mullow(C2, A2, lenA, B2, lenB, lenA + lenB - 1, ctx2));

            TIMEIT_START
            GR_MUST_SUCCEED(_gr_poly_mullow(C, A, lenA, B, lenB, lenA + lenB - 1, ctx));
            TIMEIT_STOP_VALUES(t1, __)
            TIMEIT_START
            GR_MUST_SUCCEED(_gr_poly_mullow(C2, A2, lenA, B2, lenB, lenA + lenB - 1, ctx2));
            TIMEIT_STOP_VALUES(t2, __)
            speedup_mul = t2 / t1;

            TIMEIT_START
            GR_MUST_SUCCEED(_gr_poly_divrem(Q, R, C, lenC, B, lenB, ctx));
            TIMEIT_STOP_VALUES(t1, __)
            TIMEIT_START
            GR_MUST_SUCCEED(_gr_poly_divrem(Q2, R2, C2, lenC, B2, lenB, ctx2));
            TIMEIT_STOP_VALUES(t2, __)
            speedup_divrem = t2 / t1;

            TIMEIT_START
            GR_MUST_SUCCEED(_gr_poly_mullow(C, A, lenA, A, lenA, lenA + lenB - 1, ctx));
            TIMEIT_STOP_VALUES(t1, __)
            TIMEIT_START
            GR_MUST_SUCCEED(_gr_poly_mullow(C2, A2, lenA, A2, lenA, lenA + lenB - 1, ctx2));
            TIMEIT_STOP_VALUES(t2, __)
            speedup_sqr = t2 / t1;

            TIMEIT_START
            GR_MUST_SUCCEED(_gr_poly_mullow(C, A, lenA, B, lenB, n, ctx));
            TIMEIT_STOP_VALUES(t1, __)
            TIMEIT_START
            GR_MUST_SUCCEED(_gr_poly_mullow(C2, A2, lenA, B2, lenB, n, ctx2));
            TIMEIT_STOP_VALUES(t2, __)
            speedup_mullow = t2 / t1;

            TIMEIT_START
            GR_MUST_SUCCEED(_gr_poly_mullow(C, A, lenA, A, lenA, n, ctx));
            TIMEIT_STOP_VALUES(t1, __)
            TIMEIT_START
            GR_MUST_SUCCEED(_gr_poly_mullow(C2, A2, lenA, A2, lenA, n, ctx2));
            TIMEIT_STOP_VALUES(t2, __)
            speedup_sqrlow = t2 / t1;

            if (1)
            {
                slong lenG;
                TIMEIT_START
                GR_MUST_SUCCEED(_gr_poly_gcd(Q, &lenG, A, lenA, B, lenB, ctx));
                TIMEIT_STOP_VALUES(t1, __)
                TIMEIT_START
                GR_MUST_SUCCEED(_gr_poly_gcd(Q2, &lenG, A2, lenA, B2, lenB, ctx2));
                TIMEIT_STOP_VALUES(t2, __)
                speedup_gcd = t2 / t1;
            }

            flint_printf("%5wd %5wd   %.3f  %.3f  %.3f  %.3f  %.3f  %.3f\n",
                bits, n, speedup_mul, speedup_sqr, speedup_mullow, speedup_sqrlow, speedup_divrem, speedup_gcd);

            (void) __;

            gr_heap_clear_vec(A, lenA, ctx);
            gr_heap_clear_vec(B, lenB, ctx);
            gr_heap_clear_vec(C, lenC, ctx);
            gr_heap_clear_vec(Q, lenA, ctx);
            gr_heap_clear_vec(R, lenB, ctx);

            gr_heap_clear_vec(A2, lenA, ctx2);
            gr_heap_clear_vec(B2, lenB, ctx2);
            gr_heap_clear_vec(C2, lenC, ctx2);
            gr_heap_clear_vec(Q2, lenA, ctx2);
            gr_heap_clear_vec(R2, lenB, ctx2);
        }

        gr_ctx_clear(ctx);
        gr_ctx_clear(ctx2);
    }

    fmpz_clear(p);
    flint_randclear(state);
    return 0;
}
