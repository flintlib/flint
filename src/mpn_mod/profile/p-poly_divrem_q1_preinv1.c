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

static void
mpn_mod_set_mpn2(nn_ptr res, nn_srcptr s, slong l, gr_ctx_t ctx)
{
    MPN_NORM(s, l);
    mpn_mod_set_mpn(res, s, l, ctx);
}

#define FLINT_MPN_MUL_3_2X2(R2, R1, R0, a1, a0, b1, b0) \
    do \
    { \
        ulong __tmp2, __tmp1; \
        umul_ppmm(R1, R0, a0, b0); \
        (R2) = (a1) * (b1); \
        umul_ppmm(__tmp2, __tmp1, a0, b1); \
        add_ssaaaa(R2, R1, R2, R1, __tmp2, __tmp1); \
        umul_ppmm(__tmp2, __tmp1, a1, b0); \
        add_ssaaaa(R2, R1, R2, R1, __tmp2, __tmp1); \
    } \
    while (0) \

/* original implementation to check for regressions */
int _mpn_mod_poly_divrem_q1_preinv1_old(nn_ptr Q, nn_ptr R,
                          nn_srcptr A, slong lenA, nn_srcptr B, slong lenB,
                          nn_srcptr invL, gr_ctx_t ctx)
{
    ulong q0[MPN_MOD_MAX_LIMBS];
    ulong q1[MPN_MOD_MAX_LIMBS];
    ulong t[2 * MPN_MOD_MAX_LIMBS + 1];
    ulong u[2 * MPN_MOD_MAX_LIMBS];
    slong i;
    slong nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
    int monic = mpn_mod_is_one(invL, ctx) == T_TRUE;

    if (monic)
        mpn_mod_set(q1, A + (lenA - 1) * nlimbs, ctx);
    else
        mpn_mod_mul(q1, A + (lenA - 1) * nlimbs, invL, ctx);

    mpn_mod_mul(t, q1, B + (lenB - 2) * nlimbs, ctx);
    mpn_mod_sub(t, t, A + (lenA - 2) * nlimbs, ctx);

    if (monic)
        mpn_mod_set(q0, t, ctx);
    else
        mpn_mod_mul(q0, t, invL, ctx);

    mpn_mod_mul(t, q0, B, ctx);
    mpn_mod_add(R, A, t, ctx);

    mpn_mod_neg(Q, q0, ctx);
    mpn_mod_set(Q + nlimbs, q1, ctx);
    mpn_mod_neg(q1, q1, ctx);

    if (nlimbs == 2)
    {
        slong bits = 2 * MPN_MOD_CTX_MODULUS_BITS(ctx) + 1;
        slong slimbs = (bits + FLINT_BITS - 1) / FLINT_BITS;

        if (slimbs == 3)
        {
            for (i = 1; i < lenB - 1; i++)
            {
                nn_srcptr B1ptr = B + (i - 1) * nlimbs;
                nn_srcptr Bptr = B + i * nlimbs;
                nn_srcptr Aptr = A + i * nlimbs;

                FLINT_MPN_MUL_3_2X2(t[2], t[1], t[0], q1[1], q1[0], B1ptr[1], B1ptr[0]);
                add_sssaaaaaa(t[2], t[1], t[0], t[2], t[1], t[0], 0, Aptr[1], Aptr[0]);
                FLINT_MPN_MUL_3_2X2(u[2], u[1], u[0], q0[1], q0[0], Bptr[1], Bptr[0]);
                add_sssaaaaaa(t[2], t[1], t[0], t[2], t[1], t[0], u[2], u[1], u[0]);
                mpn_mod_set_mpn2(R + i * nlimbs, t, slimbs, ctx);
            }
        }
        else
        {
            for (i = 1; i < lenB - 1; i++)
            {
                nn_srcptr B1ptr = B + (i - 1) * nlimbs;
                nn_srcptr Bptr = B + i * nlimbs;
                nn_srcptr Aptr = A + i * nlimbs;

                FLINT_MPN_MUL_2X2(t[3], t[2], t[1], t[0], q1[1], q1[0], B1ptr[1], B1ptr[0]);
                add_ssssaaaaaaaa(t[3], t[2], t[1], t[0], t[3], t[2], t[1], t[0], 0, 0, Aptr[1], Aptr[0]);
                FLINT_MPN_MUL_2X2(u[3], u[2], u[1], u[0], q0[1], q0[0], Bptr[1], Bptr[0]);
                add_sssssaaaaaaaaaa(t[4], t[3], t[2], t[1], t[0], 0, t[3], t[2], t[1], t[0], 0, u[3], u[2], u[1], u[0]);
                mpn_mod_set_mpn2(R + i * nlimbs, t, slimbs, ctx);
            }
        }
    }
    else
    {
        for (i = 1; i < lenB - 1; i++)
        {
            flint_mpn_mul_n(t, q1, B + (i - 1) * nlimbs, nlimbs);
            flint_mpn_mul_n(u, q0, B + i * nlimbs, nlimbs);
            t[2 * nlimbs] = mpn_add_n(t, t, u, 2 * nlimbs);
            ulong cy = mpn_add_n(t, t, A + i * nlimbs, nlimbs);
            mpn_add_1(t + nlimbs, t + nlimbs, nlimbs + 1, cy);
            mpn_mod_set_mpn2(R + i * nlimbs, t, 2 * nlimbs + 1, ctx);
        }
    }

    return GR_SUCCESS;
}

/* currently always slower than the algorithms in the library */
int _mpn_mod_poly_divrem_q1_preinv1_karatsuba(nn_ptr Q, nn_ptr R,
                          nn_srcptr A, slong lenA, nn_srcptr B, slong lenB,
                          nn_srcptr invL, gr_ctx_t ctx)
{
    ulong q0[MPN_MOD_MAX_LIMBS];
    ulong q1[MPN_MOD_MAX_LIMBS];
    ulong q0q1[MPN_MOD_MAX_LIMBS];
    ulong t[2 * MPN_MOD_MAX_LIMBS + 1];
    ulong u[2 * MPN_MOD_MAX_LIMBS];
    slong i;
    slong nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
    int monic = mpn_mod_is_one(invL, ctx) == T_TRUE;

    if (monic)
        mpn_mod_set(q1, A + (lenA - 1) * nlimbs, ctx);
    else
        mpn_mod_mul(q1, A + (lenA - 1) * nlimbs, invL, ctx);

    mpn_mod_mul(t, q1, B + (lenB - 2) * nlimbs, ctx);
    mpn_mod_sub(t, t, A + (lenA - 2) * nlimbs, ctx);

    if (monic)
        mpn_mod_set(q0, t, ctx);
    else
        mpn_mod_mul(q0, t, invL, ctx);

    mpn_mod_mul(t, q0, B, ctx);
    mpn_mod_add(R, A, t, ctx);

    mpn_mod_neg(Q, q0, ctx);
    mpn_mod_set(Q + nlimbs, q1, ctx);
    mpn_mod_neg(q1, q1, ctx);

    mpn_mod_add(q0q1, q0, q1, ctx);

    nn_srcptr d = MPN_MOD_CTX_MODULUS(ctx);

    for (i = 1; i < lenB - 1; i++)
    {
        if (i % 2 == 1)
        {
            flint_mpn_submod_n(R + i * nlimbs, A + i * nlimbs, t, d, nlimbs);
            flint_mpn_addmod_n(u, B + (i - 1) * nlimbs, B + i * nlimbs, d, nlimbs);
            mpn_mod_mul(t, q0q1, u, ctx);
            flint_mpn_addmod_n(R + i * nlimbs, R + i * nlimbs, t, d, nlimbs);
            mpn_mod_mul(t, q1, B + i * nlimbs, ctx);
            flint_mpn_submod_n(R + i * nlimbs, R + i * nlimbs, t, d, nlimbs);
        }
        else
        {
            flint_mpn_addmod_n(R + i * nlimbs, A + i * nlimbs, t, d, nlimbs);
            mpn_mod_mul(t, q0, B + i * nlimbs, ctx);
            flint_mpn_addmod_n(R + i * nlimbs, R + i * nlimbs, t, d, nlimbs);
        }
    }

    return GR_SUCCESS;
}


#if 1
#undef TIMEIT_END_REPEAT
#define TIMEIT_END_REPEAT(__timer, __reps) \
            } \
            timeit_stop(__timer); \
            if (__timer->cpu >= 20) \
                break; \
            __reps *= 10; \
        } \
    } while (0)
#endif

slong parameters[] = { 2, 3, 4, 6, 8, 10, 12, 14, 16, 20, 24, 32, 48, 64, 96, 128, 0 };

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

#define OLD 0
#define FMMA 1
#define FMMA_PRECOND 2
#define KARATSUBA 3
#define KARATSUBA_PRECOND 4

/*
int _mpn_mod_poly_divrem_q1_preinv1_old(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr invL, gr_ctx_t ctx);
int _mpn_mod_poly_divrem_q1_preinv1_fmma(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr invL, gr_ctx_t ctx);
int _mpn_mod_poly_divrem_q1_preinv1_fmma_precond(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr invL, gr_ctx_t ctx);
int _mpn_mod_poly_divrem_q1_preinv1_karatsuba(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr invL, gr_ctx_t ctx);
int _mpn_mod_poly_divrem_q1_preinv1_karatsuba_precond(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr invL, gr_ctx_t ctx);
int _mpn_mod_poly_divrem_q1_preinv1(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA, nn_srcptr B, slong lenB, nn_srcptr invL, gr_ctx_t ctx);
*/

#define FIND_BEST 0
#define BEST_VS_DEFAULT 1
#define OLD_VS_DEFAULT 2

void
best_table(flint_rand_t state, int comparison, gr_ctx_t ctx)
{
    gr_ptr A, B, Q, R, invL;
    double times[6], __;
    slong lenA, lenB, len, leni;
    int best, i;

    for (leni = 0; (len = parameters[leni]) != 0; leni++)
    {
        lenA = len + 1;
        lenB = len;

        A = gr_heap_init_vec(lenA, ctx);
        B = gr_heap_init_vec(lenB, ctx);
        Q = gr_heap_init_vec(2, ctx);
        R = gr_heap_init_vec(lenB - 1, ctx);
        invL = gr_heap_init_vec(1, ctx);

        randvec(A, state, lenA, ctx);
        randvec(B, state, lenB, ctx);
        GR_MUST_SUCCEED(mpn_mod_inv(invL, GR_ENTRY(B, lenB - 1, ctx->sizeof_elem), ctx));

        TIMEIT_START;
        GR_MUST_SUCCEED(_mpn_mod_poly_divrem_q1_preinv1_old(Q, R, A, lenA, B, lenB, invL, ctx));
        TIMEIT_STOP_VALUES(__, times[0]);
        TIMEIT_START;
        GR_MUST_SUCCEED(_mpn_mod_poly_divrem_q1_preinv1_fmma(Q, R, A, lenA, B, lenB, invL, ctx));
        TIMEIT_STOP_VALUES(__, times[1]);
        TIMEIT_START;
        GR_MUST_SUCCEED(_mpn_mod_poly_divrem_q1_preinv1_fmma_precond(Q, R, A, lenA, B, lenB, invL, ctx));
        TIMEIT_STOP_VALUES(__, times[2]);
        TIMEIT_START;
        GR_MUST_SUCCEED(_mpn_mod_poly_divrem_q1_preinv1_karatsuba(Q, R, A, lenA, B, lenB, invL, ctx));
        TIMEIT_STOP_VALUES(__, times[3]);
        TIMEIT_START;
        GR_MUST_SUCCEED(_mpn_mod_poly_divrem_q1_preinv1_karatsuba_precond(Q, R, A, lenA, B, lenB, invL, ctx));
        TIMEIT_STOP_VALUES(__, times[4]);
        TIMEIT_START;
        GR_MUST_SUCCEED(_mpn_mod_poly_divrem_q1_preinv1(Q, R, A, lenA, B, lenB, invL, ctx));
        TIMEIT_STOP_VALUES(__, times[5]);

        best = 0;
        for (i = 1; i < 5; i++)
        {
            if (times[i] < times[best])
                best = i;
        }

        if (comparison == FIND_BEST)
        {
            flint_printf("%6wd", best);
        }
        else if (comparison == BEST_VS_DEFAULT)
        {
            flint_printf(" %.3f", times[best] / times[5]);
        }
        else if (comparison == OLD_VS_DEFAULT)
        {
            flint_printf(" %.3f", times[0] / times[5]);
        }

        fflush(stdout);

        (void) __;

        gr_heap_clear_vec(A, lenA, ctx);
        gr_heap_clear_vec(B, lenB, ctx);
        gr_heap_clear_vec(Q, 2, ctx);
        gr_heap_clear_vec(R, lenB - 1, ctx);
    }
}

int main(void)
{
    fmpz_t p;
    gr_ctx_t ctx;
    flint_rand_t state;
    slong bits;
    slong len, leni;

    flint_rand_init(state);

    flint_printf("      ");
    for (leni = 0; (len = parameters[leni]) != 0; leni++)
        flint_printf("%5wd ", len);

    flint_printf("\n");

    for (bits = 96; bits <= 1024; bits += 32)
    {
        flint_printf("%5wd", bits);
        fflush(stdout);

        fmpz_init(p);
        fmpz_randprime(p, state, bits, 0);
        GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(ctx, p));
        //best_table(state, BEST_VS_DEFAULT, ctx);
        //best_table(state, FIND_BEST, ctx);
        best_table(state, OLD_VS_DEFAULT, ctx);
        gr_ctx_clear(ctx);

        flint_printf("\n");
    }

    fmpz_clear(p);
    flint_rand_clear(state);
}

