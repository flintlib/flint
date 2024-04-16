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
            if (__timer->cpu >= 20) \
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

#define MUL 0
#define SQR 1
#define MULLOW 2
#define SQRLOW 3
#define UNBALANCED 4

slong cutoff_bits[1000];
slong karatsuba_cutoffs[1000];
slong KS_or_fft_small_cutoffs[1000];
slong num_cutoffs = 0;

void
cutoffs(int operation, flint_rand_t state, gr_ctx_t ctx)
{
    gr_ptr A, B, C;
    double t_classical, t_karatsuba, t_KS, t_fft_small, best, __;
    slong lenA, lenB, len, n, ni;
    int squaring;
    char * beststr = "none";
    slong bits = MPN_MOD_CTX_MODULUS_BITS(ctx);

    int found_nonclassical = 0;

    cutoff_bits[num_cutoffs] = bits;

    for (n = 4; ; n++)
    {
        lenA = n;
        lenB = (operation == UNBALANCED) ? 8 * n : n;
        len = (operation == MULLOW || operation == SQRLOW) ? n : lenA + lenB - 1;
        squaring = (operation == SQR || operation == SQRLOW);

        A = gr_heap_init_vec(lenA, ctx);
        B = gr_heap_init_vec(lenB, ctx);
        C = gr_heap_init_vec(len, ctx);

        randvec(A, state, lenA, ctx);
        randvec(B, state, lenB, ctx);

        TIMEIT_START
        GR_MUST_SUCCEED(_mpn_mod_poly_mullow_classical(C, A, lenA, squaring ? A : B, lenB, len, ctx));
        TIMEIT_STOP_VALUES(t_classical, __)

        TIMEIT_START
        GR_MUST_SUCCEED(_mpn_mod_poly_mullow_karatsuba(C, A, lenA, squaring ? A : B, lenB, len, -1, ctx));
        TIMEIT_STOP_VALUES(t_karatsuba, __)

        TIMEIT_START
        GR_MUST_SUCCEED(_mpn_mod_poly_mullow_KS(C, A, lenA, squaring ? A : B, lenB, len, ctx));
        TIMEIT_STOP_VALUES(t_KS, __)

        if (_mpn_mod_poly_mullow_fft_small(C, A, lenA, squaring ? A : B, lenB, len, ctx) == GR_SUCCESS)
        {
            TIMEIT_START
            GR_MUST_SUCCEED(_mpn_mod_poly_mullow_fft_small(C, A, lenA, squaring ? A : B, lenB, len, ctx));
            TIMEIT_STOP_VALUES(t_fft_small, __)
        }
        else
        {
            t_fft_small = D_INF;
        }

        best = FLINT_MIN(t_classical, t_karatsuba);
        best = FLINT_MIN(best, t_KS);
        best = FLINT_MIN(best, t_fft_small);

        if (best == t_fft_small) beststr = "fft_small";
        if (best == t_KS) beststr = "KS";
        if (best == t_karatsuba) beststr = "karatsuba";
        if (best == t_classical) beststr = "classical";

        flint_printf("%wd   %wd    %g  %g  %g  %g     %s\n", bits, n, t_classical, t_karatsuba, t_KS, t_fft_small, beststr);

        (void) __;

        gr_heap_clear_vec(A, lenA, ctx);
        gr_heap_clear_vec(B, lenB, ctx);
        gr_heap_clear_vec(C, len, ctx);

        if (!found_nonclassical)
        {
            if (best != t_classical)
            {
                karatsuba_cutoffs[num_cutoffs] = n;
                found_nonclassical = 1;
            }
        }

        if (best != t_classical && best != t_karatsuba && (t_fft_small == D_INF || (t_fft_small < t_classical && t_fft_small < t_karatsuba)))
        {
            KS_or_fft_small_cutoffs[num_cutoffs] = n;

            num_cutoffs++;
            break;
        }
    }

    for (ni = 0; ni < num_cutoffs; ni++)
        flint_printf("  {%wd, %wd},   /* bits = %wd */\n", karatsuba_cutoffs[ni], KS_or_fft_small_cutoffs[ni], cutoff_bits[ni]);
}

void
compare(int operation, slong nmax, flint_rand_t state, gr_ctx_t ctx)
{
    gr_ptr A, B, C;
    double t_classical, t_karatsuba, t_KS, t_fft_small, t_default, best, __;
    slong lenA, lenB, len;
    int squaring;
    char * beststr = "none";
    slong bits = MPN_MOD_CTX_MODULUS_BITS(ctx);

    lenA = 1 + n_randint(state, 1 + n_randint(state, nmax));
    lenB = (n_randint(state, 2) || operation == SQR || operation == SQRLOW) ? lenA : 1 + n_randint(state, 1 + n_randint(state, nmax));
    len = (n_randint(state, 2) || (operation != MULLOW && operation != SQRLOW)) ? lenA + lenB - 1 : 1 + n_randint(state, lenA + lenB - 1);
    squaring = (operation == SQR || operation == SQRLOW);

    A = gr_heap_init_vec(lenA, ctx);
    B = gr_heap_init_vec(lenB, ctx);
    C = gr_heap_init_vec(len, ctx);

    randvec(A, state, lenA, ctx);
    randvec(B, state, lenB, ctx);

    TIMEIT_START
    GR_MUST_SUCCEED(_mpn_mod_poly_mullow(C, A, lenA, squaring ? A : B, lenB, len, ctx));
    TIMEIT_STOP_VALUES(t_default, __)

    TIMEIT_START
    GR_MUST_SUCCEED(_mpn_mod_poly_mullow_classical(C, A, lenA, squaring ? A : B, lenB, len, ctx));
    TIMEIT_STOP_VALUES(t_classical, __)

    TIMEIT_START
    GR_MUST_SUCCEED(_mpn_mod_poly_mullow_karatsuba(C, A, lenA, squaring ? A : B, lenB, len, -1, ctx));
    TIMEIT_STOP_VALUES(t_karatsuba, __)

    TIMEIT_START
    GR_MUST_SUCCEED(_mpn_mod_poly_mullow_KS(C, A, lenA, squaring ? A : B, lenB, len, ctx));
    TIMEIT_STOP_VALUES(t_KS, __)

    if (_mpn_mod_poly_mullow_fft_small(C, A, lenA, squaring ? A : B, lenB, len, ctx) == GR_SUCCESS)
    {
        TIMEIT_START
        GR_MUST_SUCCEED(_mpn_mod_poly_mullow_fft_small(C, A, lenA, squaring ? A : B, lenB, len, ctx));
        TIMEIT_STOP_VALUES(t_fft_small, __)
    }
    else
    {
        t_fft_small = D_INF;
    }

    best = FLINT_MIN(t_classical, t_karatsuba);
    best = FLINT_MIN(best, t_KS);
    best = FLINT_MIN(best, t_fft_small);

    if (best == t_fft_small) beststr = "fft_small";
    if (best == t_KS) beststr = "KS";
    if (best == t_karatsuba) beststr = "karatsuba";
    if (best == t_classical) beststr = "classical";

    flint_printf("%5wd    %5wd %5wd %5wd    %-1.6f  %-1.6f  %-1.6f  %-1.6f     %10s    %-1.6f    %.3f\n",
        bits, lenA, lenB, len, t_classical, t_karatsuba, t_KS, t_fft_small == D_INF ? 0.0 : t_fft_small, beststr, t_default, t_default / best);

    (void) __;

    gr_heap_clear_vec(A, lenA, ctx);
    gr_heap_clear_vec(B, lenB, ctx);
    gr_heap_clear_vec(C, len, ctx);

}

int main()
{
    fmpz_t p;
    gr_ctx_t ctx;
    flint_rand_t state;
    slong bits;

    flint_randinit(state);

    for (bits = 1; bits <= 100; bits++)
    {
        gr_ctx_init_mpn_mod_randtest(ctx, state);
        compare(MUL, 200, state, ctx);
        gr_ctx_clear(ctx);
    }

    for (bits = 80; bits <= 1024; bits += 16)
    {
        fmpz_init(p);
        fmpz_randprime(p, state, bits, 0);
        GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(ctx, p));
        cutoffs(MUL, state, ctx);
        gr_ctx_clear(ctx);
    }

    fmpz_clear(p);
    flint_randclear(state);
}
