/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_mod.h"

static int _mpn_mod_poly_divrem_q1_preinv1_fmma_n(nn_ptr Q, nn_ptr R,
                          nn_srcptr A, slong lenA, nn_srcptr B, slong lenB,
                          nn_srcptr invL, gr_ctx_t ctx)
{
    ulong q0[MPN_MOD_MAX_LIMBS];
    ulong q1[MPN_MOD_MAX_LIMBS];
    ulong t[MPN_MOD_MAX_LIMBS];
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

    nn_srcptr d = MPN_MOD_CTX_MODULUS(ctx);
    flint_bitcnt_t norm = MPN_MOD_CTX_NORM(ctx);
    nn_srcptr dnormed = MPN_MOD_CTX_MODULUS_NORMED(ctx);
    nn_srcptr dinv = MPN_MOD_CTX_MODULUS_PREINV(ctx);

    for (i = 1; i < lenB - 1; i++)
    {
        flint_mpn_fmmamod_preinvn(t, q1, B + (i - 1) * nlimbs, q0, B + i * nlimbs, nlimbs, dnormed, dinv, norm);
        flint_mpn_addmod_n(R + i * nlimbs, A + i * nlimbs, t, d, nlimbs);
    }

    return GR_SUCCESS;
}

int _mpn_mod_poly_divrem_q1_preinv1_fmma(nn_ptr Q, nn_ptr R,
                          nn_srcptr A, slong lenA, nn_srcptr B, slong lenB,
                          nn_srcptr invL, gr_ctx_t ctx)
{
    ulong q0[2];
    ulong q1[2];
    ulong t[2];
    slong i;
    slong nlimbs = MPN_MOD_CTX_NLIMBS(ctx);

    if (nlimbs != 2)
        return _mpn_mod_poly_divrem_q1_preinv1_fmma_n(Q, R, A, lenA, B, lenB, invL, ctx);

    int monic = (invL[0] == 1 && invL[1] == 0);

    nn_srcptr d = MPN_MOD_CTX_MODULUS(ctx);
    flint_bitcnt_t norm = MPN_MOD_CTX_NORM(ctx);
    nn_srcptr dnormed = MPN_MOD_CTX_MODULUS_NORMED(ctx);
    nn_srcptr dinv = MPN_MOD_CTX_MODULUS_PREINV(ctx);

    if (monic)
        flint_mpn_copyi(q1, A + (lenA - 1) * nlimbs, 2);
    else
        flint_mpn_mulmod_preinvn_2(q1, A + (lenA - 1) * nlimbs, invL, dnormed, dinv, norm);

    flint_mpn_mulmod_preinvn_2(t, q1, B + (lenB - 2) * nlimbs, dnormed, dinv, norm);
    flint_mpn_submod_2(t, t, A + (lenA - 2) * nlimbs, d);

    if (monic)
        flint_mpn_copyi(q0, t, 2);
    else
        flint_mpn_mulmod_preinvn_2(q0, t, invL, dnormed, dinv, norm);

    flint_mpn_mulmod_preinvn_2(t, q0, B, dnormed, dinv, norm);
    flint_mpn_addmod_2(R, A, t, d);
    flint_mpn_negmod_2(Q, q0, d);
    flint_mpn_copyi(Q + nlimbs, q1, 2);
    flint_mpn_negmod_2(q1, q1, d);

    for (i = 1; i < lenB - 1; i++)
    {
        flint_mpn_fmmamod_preinvn_2(t, q1, B + (i - 1) * nlimbs, q0, B + i * nlimbs, dnormed, dinv, norm);
        flint_mpn_addmod_2(R + i * nlimbs, A + i * nlimbs, t, d);
    }

    return GR_SUCCESS;
}

int _mpn_mod_poly_divrem_q1_preinv1_fmma_precond(nn_ptr Q, nn_ptr R,
                          nn_srcptr A, slong lenA, nn_srcptr B, slong lenB,
                          nn_srcptr invL, gr_ctx_t ctx)
{
    ulong q0[MPN_MOD_MAX_LIMBS];
    ulong q1[MPN_MOD_MAX_LIMBS];
    ulong t[2 * MPN_MOD_MAX_LIMBS + 1];
    slong i;
    slong nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
    int monic = mpn_mod_is_one(invL, ctx) == T_TRUE;

    ulong tq0[MPN_MOD_MAX_LIMBS * MPN_MOD_MAX_LIMBS + 1];
    ulong tq1[MPN_MOD_MAX_LIMBS * MPN_MOD_MAX_LIMBS + 1];

    if (monic)
        mpn_mod_set(q1, A + (lenA - 1) * nlimbs, ctx);
    else
        mpn_mod_mul(q1, A + (lenA - 1) * nlimbs, invL, ctx);

    mpn_mod_set(Q + nlimbs, q1, ctx);
    mpn_mod_neg(q1, q1, ctx);

    nn_srcptr d = MPN_MOD_CTX_MODULUS(ctx);
    flint_bitcnt_t norm = MPN_MOD_CTX_NORM(ctx);
    nn_srcptr dnormed = MPN_MOD_CTX_MODULUS_NORMED(ctx);
    nn_srcptr dinv = MPN_MOD_CTX_MODULUS_PREINV(ctx);

    flint_mpn_mulmod_precond_matrix_precompute(tq1, q1, nlimbs, dnormed, dinv, norm);
    flint_mpn_mulmod_precond_matrix(t, tq1, B + (lenB - 2) * nlimbs, nlimbs, dnormed, dinv, norm);
    mpn_mod_add(t, t, A + (lenA - 2) * nlimbs, ctx);

    if (monic)
    {
        mpn_mod_neg(q0, t, ctx);
    }
    else
    {
        mpn_mod_mul(q0, t, invL, ctx);
        mpn_mod_neg(q0, q0, ctx);
    }

    flint_mpn_mulmod_precond_matrix_precompute(tq0, q0, nlimbs, dnormed, dinv, norm);
    flint_mpn_mulmod_precond_matrix(t, tq0, B, nlimbs, dnormed, dinv, norm);
    mpn_mod_add(R, A, t, ctx);
    mpn_mod_neg(Q, q0, ctx);

    for (i = 1; i < lenB - 1; i++)
    {
        flint_mpn_fmmamod_precond_matrix(t, tq1, B + (i - 1) * nlimbs, tq0, B + i * nlimbs, nlimbs, dnormed, dinv, norm);
        flint_mpn_addmod_n(R + i * nlimbs, A + i * nlimbs, t, d, nlimbs);
    }

    return GR_SUCCESS;
}

int _mpn_mod_poly_divrem_q1_preinv1_karatsuba_precond(nn_ptr Q, nn_ptr R,
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

    ulong tq0[MPN_MOD_MAX_LIMBS * MPN_MOD_MAX_LIMBS + 1];
    ulong tq1[MPN_MOD_MAX_LIMBS * MPN_MOD_MAX_LIMBS + 1];
    ulong tq0q1[MPN_MOD_MAX_LIMBS * MPN_MOD_MAX_LIMBS + 1];

    if (monic)
        mpn_mod_set(q1, A + (lenA - 1) * nlimbs, ctx);
    else
        mpn_mod_mul(q1, A + (lenA - 1) * nlimbs, invL, ctx);

    mpn_mod_set(Q + nlimbs, q1, ctx);
    mpn_mod_neg(q1, q1, ctx);

    flint_bitcnt_t norm = MPN_MOD_CTX_NORM(ctx);
    nn_srcptr dnormed = MPN_MOD_CTX_MODULUS_NORMED(ctx);
    nn_srcptr dinv = MPN_MOD_CTX_MODULUS_PREINV(ctx);

    flint_mpn_mulmod_precond_matrix_precompute(tq1, q1, nlimbs, dnormed, dinv, norm);
    flint_mpn_mulmod_precond_matrix(t, tq1, B + (lenB - 2) * nlimbs, nlimbs, dnormed, dinv, norm);

    mpn_mod_add(t, t, A + (lenA - 2) * nlimbs, ctx);

    if (monic)
    {
        mpn_mod_neg(q0, t, ctx);
    }
    else
    {
        mpn_mod_mul(q0, t, invL, ctx);
        mpn_mod_neg(q0, q0, ctx);
    }

    /* mpn_mod_mul(t, q0, B, ctx); */
    flint_mpn_mulmod_precond_matrix_precompute(tq0, q0, nlimbs, dnormed, dinv, norm);
    flint_mpn_mulmod_precond_matrix(t, tq0, B, nlimbs, dnormed, dinv, norm);

    mpn_mod_add(R, A, t, ctx);
    mpn_mod_neg(Q, q0, ctx);
    mpn_mod_add(q0q1, q0, q1, ctx);

    /* flint_mpn_mulmod_precond_precompute(tq0q1, q0q1, nlimbs, dnormed, dinv, norm); */
    for (i = 0; i < nlimbs; i++)
        flint_mpn_addmod_n(tq0q1 + i * nlimbs, tq0 + i * nlimbs, tq1 + i * nlimbs, dnormed, nlimbs);

    nn_srcptr d = MPN_MOD_CTX_MODULUS(ctx);

    for (i = 1; i < lenB - 1; i++)
    {
        if (i % 2 == 1)
        {
            flint_mpn_submod_n(R + i * nlimbs, A + i * nlimbs, t, d, nlimbs);
            flint_mpn_addmod_n(u, B + (i - 1) * nlimbs, B + i * nlimbs, d, nlimbs);
            flint_mpn_mulmod_precond_matrix(t, tq0q1, u, nlimbs, dnormed, dinv, norm);
            flint_mpn_addmod_n(R + i * nlimbs, R + i * nlimbs, t, d, nlimbs);
            flint_mpn_mulmod_precond_matrix(t, tq1, B + i * nlimbs, nlimbs, dnormed, dinv, norm);
            flint_mpn_submod_n(R + i * nlimbs, R + i * nlimbs, t, d, nlimbs);
        }
        else
        {
            flint_mpn_addmod_n(R + i * nlimbs, A + i * nlimbs, t, d, nlimbs);
            flint_mpn_mulmod_precond_matrix(t, tq0, B + i * nlimbs, nlimbs, dnormed, dinv, norm);
            flint_mpn_addmod_n(R + i * nlimbs, R + i * nlimbs, t, d, nlimbs);
        }
    }

    return GR_SUCCESS;
}

/* tuned using p-poly_divrem_q1_preinv1 */
int _mpn_mod_poly_divrem_q1_preinv1(nn_ptr Q, nn_ptr R,
                          nn_srcptr A, slong lenA, nn_srcptr B, slong lenB,
                          nn_srcptr invL, gr_ctx_t ctx)
{
    slong bits = MPN_MOD_CTX_MODULUS_BITS(ctx);

    FLINT_ASSERT(lenB >= 2);

    if (lenB < 10 || bits <= 576 || (bits <= 704 && lenB < 14))
        return _mpn_mod_poly_divrem_q1_preinv1_fmma(Q, R, A, lenA, B, lenB, invL, ctx);

    if (bits > 896 || (bits > 768 && lenB > 12) || lenB > 20)
        return _mpn_mod_poly_divrem_q1_preinv1_karatsuba_precond(Q, R, A, lenA, B, lenB, invL, ctx);
    else
        return _mpn_mod_poly_divrem_q1_preinv1_fmma_precond(Q, R, A, lenA, B, lenB, invL, ctx);
}

