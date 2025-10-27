/*
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_mod.h"

int
_mpn_mod_vec_clear(nn_ptr FLINT_UNUSED(res), slong FLINT_UNUSED(len), gr_ctx_t FLINT_UNUSED(ctx))
{
    return GR_SUCCESS;
}

int
_mpn_mod_vec_zero(nn_ptr res, slong len, gr_ctx_t ctx)
{
    flint_mpn_zero(res, len * MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

int
_mpn_mod_vec_set(nn_ptr res, nn_srcptr x, slong len, gr_ctx_t ctx)
{
    flint_mpn_copyi(res, x, len * MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

void
_mpn_mod_vec_swap(nn_ptr vec1, nn_ptr vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    slong n = MPN_MOD_CTX_NLIMBS(ctx);
    for (i = 0; i < len * n; i++)
        FLINT_SWAP(ulong, vec1[i], vec2[i]);
}

int
_mpn_mod_vec_neg(nn_ptr res, nn_srcptr x, slong len, gr_ctx_t ctx)
{
    slong n = MPN_MOD_CTX_NLIMBS(ctx);
    nn_srcptr d = MPN_MOD_CTX_MODULUS(ctx);
    slong i;

    if (n == 2)
    {
        /* Only read to registers once */
        ulong dd[2];
        dd[0] = d[0];
        dd[1] = d[1];

        for (i = 0; i < len; i++)
            flint_mpn_negmod_2(res + i * n, x + i * n, dd);
    }
    else
        for (i = 0; i < len; i++)
            flint_mpn_negmod_n(res + i * n, x + i * n, d, n);

    return GR_SUCCESS;
}

int
_mpn_mod_vec_add(nn_ptr res, nn_srcptr x, nn_srcptr y, slong len, gr_ctx_t ctx)
{
    slong n = MPN_MOD_CTX_NLIMBS(ctx);
    nn_srcptr d = MPN_MOD_CTX_MODULUS(ctx);
    slong i;

    if (n == 2)
    {
        /* Only read to registers once */
        ulong dd[2];
        dd[0] = d[0];
        dd[1] = d[1];

        if (MPN_MOD_CTX_NORM(ctx))
            for (i = 0; i < len; i++)
                _flint_mpn_addmod_2(res + i * n, x + i * n, y + i * n, dd);
        else
            for (i = 0; i < len; i++)
                flint_mpn_addmod_2(res + i * n, x + i * n, y + i * n, dd);
    }
    else
        for (i = 0; i < len; i++)
            flint_mpn_addmod_n(res + i * n, x + i * n, y + i * n, d, n);

    return GR_SUCCESS;
}

int
_mpn_mod_vec_sub(nn_ptr res, nn_srcptr x, nn_srcptr y, slong len, gr_ctx_t ctx)
{
    slong n = MPN_MOD_CTX_NLIMBS(ctx);
    nn_srcptr d = MPN_MOD_CTX_MODULUS(ctx);
    slong i;

    if (n == 2)
    {
        /* Only read to registers once */
        ulong dd[2];
        dd[0] = d[0];
        dd[1] = d[1];

        for (i = 0; i < len; i++)
            flint_mpn_submod_2(res + i * n, x + i * n, y + i * n, dd);
    }
    else
        for (i = 0; i < len; i++)
            flint_mpn_submod_n(res + i * n, x + i * n, y + i * n, d, n);

    return GR_SUCCESS;
}

int
_mpn_mod_vec_mul(nn_ptr res, nn_srcptr x, nn_srcptr y, slong len, gr_ctx_t ctx)
{
    slong n = MPN_MOD_CTX_NLIMBS(ctx);
    slong i;

    if (n == 2)
    {
        flint_bitcnt_t norm = MPN_MOD_CTX_NORM(ctx);
        nn_srcptr dnormed = MPN_MOD_CTX_MODULUS_NORMED(ctx);
        nn_srcptr dinv = MPN_MOD_CTX_MODULUS_PREINV(ctx);

        for (i = 0; i < len; i++)
            flint_mpn_mulmod_preinvn_2(res + i * n, x + i * n, y + i * n, dnormed, dinv, norm);
    }
    else
        for (i = 0; i < len; i++)
            mpn_mod_mul(res + i * n, x + i * n, y + i * n, ctx);

    return GR_SUCCESS;
}

static void
flint_mpn_mulmod_precond_shoup_2_norm1(mp_ptr res, ulong a1, ulong a0, ulong apre1, ulong apre0, mp_srcptr b, ulong d1, ulong d0)
{
    mp_limb_t t3, t2, t1, t0;

    FLINT_MPN_MUL_2X2(t3, t2, t1, t0, apre1, apre0, b[1], b[0]);
    FLINT_MPN_MULLOW_2X2(t1, t0, t3, t2, d1, d0);
    FLINT_MPN_MULLOW_2X2(t3, t2, a1, a0, b[1], b[0]);
    sub_ddmmss(t3, t2, t3, t2, t1, t0);
    if (t3 > d1 || (t3 == d1 && t2 >= d0))
        sub_ddmmss(t3, t2, t3, t2, d1, d0);
    res[0] = t2;
    res[1] = t3;
}

static void
flint_mpn_mulmod_precond_shoup_2_norm0(mp_ptr res, ulong a1, ulong a0, ulong apre1, ulong apre0, mp_srcptr b, ulong d1, ulong d0)
{
    mp_limb_t cy, tcy, ucy, t3, t2, t1, t0;

    FLINT_MPN_MUL_2X2(t3, t2, t1, t0, apre1, apre0, b[1], b[0]);
    FLINT_MPN_MUL_3P2X2(ucy, t1, t0, t3, t2, d1, d0);
    FLINT_MPN_MUL_3P2X2(tcy, t3, t2, a1, a0, b[1], b[0]);
    sub_dddmmmsss(cy, t3, t2, tcy, t3, t2, ucy, t1, t0);
    if (cy || (t3 > d1 || (t3 == d1 && t2 >= d0)))
        sub_ddmmss(t3, t2, t3, t2, d1, d0);
    res[0] = t2;
    res[1] = t3;
}

static int
want_shoup_2(slong len, ulong norm)
{
    return len >= (norm ? 4 : 6);
}

/* todo: worth it to check for special cases (0, 1)? */
int
_mpn_mod_vec_mul_scalar(nn_ptr res, nn_srcptr x, slong len, nn_srcptr y, gr_ctx_t ctx)
{
    slong n = MPN_MOD_CTX_NLIMBS(ctx);
    slong i;
    int precond;

    flint_bitcnt_t norm = MPN_MOD_CTX_NORM(ctx);
    nn_srcptr dnormed = MPN_MOD_CTX_MODULUS_NORMED(ctx);
    nn_srcptr dinv = MPN_MOD_CTX_MODULUS_PREINV(ctx);
    nn_srcptr d = MPN_MOD_CTX_MODULUS(ctx);

    if (n == 2)
    {
        if (!want_shoup_2(len, norm))
        {
            for (i = 0; i < len; i++)
                flint_mpn_mulmod_preinvn_2(res + i * n, x + i * n, y, dnormed, dinv, norm);
        }
        else
        {
            ulong ypre[2];
            ulong y0, y1, ypre0, ypre1, d0, d1;

            flint_mpn_mulmod_precond_shoup_precompute(ypre, y, n, dnormed, dinv, norm);

            y0 = y[0];
            y1 = y[1];
            ypre0 = ypre[0];
            ypre1 = ypre[1];
            d0 = d[0];
            d1 = d[1];

            if (norm == 0)
                for (i = 0; i < len; i++)
                    flint_mpn_mulmod_precond_shoup_2_norm0(res + i * n, y1, y0, ypre1, ypre0, x + i * n, d1, d0);
            else
                for (i = 0; i < len; i++)
                    flint_mpn_mulmod_precond_shoup_2_norm1(res + i * n, y1, y0, ypre1, ypre0, x + i * n, d1, d0);
        }
    }
    else
    {
        precond = flint_mpn_mulmod_want_precond(n, len, norm);

        if (precond == MPN_MULMOD_PRECOND_SHOUP)
        {
            TMP_INIT;
            TMP_START;
            nn_ptr ypre = TMP_ALLOC(n * sizeof(mp_limb_t));
            flint_mpn_mulmod_precond_shoup_precompute(ypre, y, n, dnormed, dinv, norm);
            for (i = 0; i < len; i++)
                flint_mpn_mulmod_precond_shoup(res + i * n, y, ypre, x + i * n, n, d, norm);
            TMP_END;
        }
        else if (precond == MPN_MULMOD_PRECOND_MATRIX)
        {
            TMP_INIT;
            TMP_START;
            nn_ptr ypre = TMP_ALLOC(flint_mpn_mulmod_precond_matrix_alloc(n) * sizeof(mp_limb_t));
            flint_mpn_mulmod_precond_matrix_precompute(ypre, y, n, dnormed, dinv, norm);
            for (i = 0; i < len; i++)
                flint_mpn_mulmod_precond_matrix(res + i * n, ypre, x + i * n, n, dnormed, dinv, norm);
            TMP_END;
        }
        else
        {
            for (i = 0; i < len; i++)
                mpn_mod_mul(res + i * n, x + i * n, y, ctx);
        }
    }

    return GR_SUCCESS;
}

int
_mpn_mod_scalar_mul_vec(nn_ptr res, nn_srcptr y, nn_srcptr x, slong len, gr_ctx_t ctx)
{
    return _mpn_mod_vec_mul_scalar(res, x, len, y, ctx);
}

/* todo: worth it to check for special cases (0, 1)? */
/* todo: shoup multiplication */
int
_mpn_mod_vec_addmul_scalar(nn_ptr res, nn_srcptr x, slong len, nn_srcptr y, gr_ctx_t ctx)
{
    slong n = MPN_MOD_CTX_NLIMBS(ctx);
    slong i;
    int precond;

    flint_bitcnt_t norm = MPN_MOD_CTX_NORM(ctx);
    nn_srcptr dnormed = MPN_MOD_CTX_MODULUS_NORMED(ctx);
    nn_srcptr dinv = MPN_MOD_CTX_MODULUS_PREINV(ctx);
    nn_srcptr d = MPN_MOD_CTX_MODULUS(ctx);

    if (n == 2)
    {
        if (!want_shoup_2(len, norm))
        {
            ulong t[2];

            for (i = 0; i < len; i++)
            {
                flint_mpn_mulmod_preinvn_2(t, x + i * n, y, dnormed, dinv, norm);
                flint_mpn_addmod_2(res + i * n, res + i * n, t, d);
            }
        }
        else
        {
            ulong ypre[2];
            ulong y0, y1, ypre0, ypre1, d0, d1;
            ulong t[2];

            flint_mpn_mulmod_precond_shoup_precompute(ypre, y, n, dnormed, dinv, norm);

            y0 = y[0];
            y1 = y[1];
            ypre0 = ypre[0];
            ypre1 = ypre[1];
            d0 = d[0];
            d1 = d[1];

            if (norm == 0)
            {
                for (i = 0; i < len; i++)
                {
                    flint_mpn_mulmod_precond_shoup_2_norm0(t, y1, y0, ypre1, ypre0, x + i * n, d1, d0);
                    flint_mpn_addmod_2(res + i * n, res + i * n, t, d);
                }
            }
            else
            {
                for (i = 0; i < len; i++)
                {
                    flint_mpn_mulmod_precond_shoup_2_norm1(t, y1, y0, ypre1, ypre0, x + i * n, d1, d0);
                    flint_mpn_addmod_2(res + i * n, res + i * n, t, d);
                }
            }
        }
    }
    else
    {
        ulong t[MPN_MOD_MAX_LIMBS];

        precond = flint_mpn_mulmod_want_precond(n, len, norm);

        if (precond == MPN_MULMOD_PRECOND_SHOUP)
        {
            TMP_INIT;
            TMP_START;
            nn_ptr ypre = TMP_ALLOC(n * sizeof(mp_limb_t));
            flint_mpn_mulmod_precond_shoup_precompute(ypre, y, n, dnormed, dinv, norm);
            for (i = 0; i < len; i++)
            {
                flint_mpn_mulmod_precond_shoup(t, y, ypre, x + i * n, n, d, norm);
                mpn_mod_add(res + i * n, res + i * n, t, ctx);
            }
            TMP_END;
        }
        else if (precond == MPN_MULMOD_PRECOND_MATRIX)
        {
            TMP_INIT;
            TMP_START;
            nn_ptr ypre = TMP_ALLOC(flint_mpn_mulmod_precond_matrix_alloc(n) * sizeof(mp_limb_t));
            flint_mpn_mulmod_precond_matrix_precompute(ypre, y, n, dnormed, dinv, norm);
            for (i = 0; i < len; i++)
            {
                flint_mpn_mulmod_precond_matrix(t, ypre, x + i * n, n, dnormed, dinv, norm);
                mpn_mod_add(res + i * n, res + i * n, t, ctx);
            }
            TMP_END;
        }
        else
        {
            for (i = 0; i < len; i++)
            {
                mpn_mod_mul(t, x + i * n, y, ctx);
                mpn_mod_add(res + i * n, res + i * n, t, ctx);
            }
        }
    }

    return GR_SUCCESS;
}

int
_mpn_mod_vec_submul_scalar(nn_ptr res, nn_srcptr x, slong len, nn_srcptr y, gr_ctx_t ctx)
{
    slong n = MPN_MOD_CTX_NLIMBS(ctx);
    slong i;
    int precond;

    flint_bitcnt_t norm = MPN_MOD_CTX_NORM(ctx);
    nn_srcptr dnormed = MPN_MOD_CTX_MODULUS_NORMED(ctx);
    nn_srcptr dinv = MPN_MOD_CTX_MODULUS_PREINV(ctx);
    nn_srcptr d = MPN_MOD_CTX_MODULUS(ctx);

    if (n == 2)
    {
        if (!want_shoup_2(len, norm))
        {
            ulong t[2];

            for (i = 0; i < len; i++)
            {
                flint_mpn_mulmod_preinvn_2(t, x + i * n, y, dnormed, dinv, norm);
                flint_mpn_submod_2(res + i * n, res + i * n, t, d);
            }
        }
        else
        {
            ulong ypre[2];
            ulong y0, y1, ypre0, ypre1, d0, d1;
            ulong t[2];

            flint_mpn_mulmod_precond_shoup_precompute(ypre, y, n, dnormed, dinv, norm);

            y0 = y[0];
            y1 = y[1];
            ypre0 = ypre[0];
            ypre1 = ypre[1];
            d0 = d[0];
            d1 = d[1];

            if (norm == 0)
            {
                for (i = 0; i < len; i++)
                {
                    flint_mpn_mulmod_precond_shoup_2_norm0(t, y1, y0, ypre1, ypre0, x + i * n, d1, d0);
                    flint_mpn_submod_2(res + i * n, res + i * n, t, d);
                }
            }
            else
            {
                for (i = 0; i < len; i++)
                {
                    flint_mpn_mulmod_precond_shoup_2_norm1(t, y1, y0, ypre1, ypre0, x + i * n, d1, d0);
                    flint_mpn_submod_2(res + i * n, res + i * n, t, d);
                }
            }
        }
    }
    else
    {
        ulong t[MPN_MOD_MAX_LIMBS];

        precond = flint_mpn_mulmod_want_precond(n, len, norm);

        if (precond == MPN_MULMOD_PRECOND_SHOUP)
        {
            TMP_INIT;
            TMP_START;
            nn_ptr ypre = TMP_ALLOC(n * sizeof(mp_limb_t));
            flint_mpn_mulmod_precond_shoup_precompute(ypre, y, n, dnormed, dinv, norm);
            for (i = 0; i < len; i++)
            {
                flint_mpn_mulmod_precond_shoup(t, y, ypre, x + i * n, n, d, norm);
                mpn_mod_sub(res + i * n, res + i * n, t, ctx);
            }
            TMP_END;
        }
        else if (precond == MPN_MULMOD_PRECOND_MATRIX)
        {
            TMP_INIT;
            TMP_START;
            nn_ptr ypre = TMP_ALLOC(flint_mpn_mulmod_precond_matrix_alloc(n) * sizeof(mp_limb_t));
            flint_mpn_mulmod_precond_matrix_precompute(ypre, y, n, dnormed, dinv, norm);
            for (i = 0; i < len; i++)
            {
                flint_mpn_mulmod_precond_matrix(t, ypre, x + i * n, n, dnormed, dinv, norm);
                mpn_mod_sub(res + i * n, res + i * n, t, ctx);
            }
            TMP_END;
        }
        else
        {
            for (i = 0; i < len; i++)
            {
                mpn_mod_mul(t, x + i * n, y, ctx);
                mpn_mod_sub(res + i * n, res + i * n, t, ctx);
            }
        }
    }

    return GR_SUCCESS;
}

#include "gr.h"

/* todo: optimize for when 2n rather than 2n+1 limbs suffice */
int
_mpn_mod_vec_dot(nn_ptr res, nn_srcptr initial, int subtract, nn_srcptr vec1, nn_srcptr vec2, slong len, gr_ctx_t ctx)
{
    ulong s[2 * MPN_MOD_MAX_LIMBS + 1];
    ulong t[2 * MPN_MOD_MAX_LIMBS];
    slong n = MPN_MOD_CTX_NLIMBS(ctx);
    slong sn;
    slong i;

    if (len <= 2)
    {
        if (len == 0)
        {
            if (initial == NULL)
                flint_mpn_zero(res, n);
            else
                flint_mpn_copyi(res, initial, n);
        }
        else
        {
            if (initial == NULL)
            {
                if (len == 1)
                    mpn_mod_mul(res, vec1, vec2, ctx);
                else
                    mpn_mod_fmma(res, vec1, vec2, vec1 + n, vec2 + n, ctx);

                if (subtract)
                    mpn_mod_neg(res, res, ctx);
            }
            else
            {
                if (len == 1)
                    mpn_mod_mul(t, vec1, vec2, ctx);
                else
                    mpn_mod_fmma(t, vec1, vec2, vec1 + n, vec2 + n, ctx);

                if (subtract)
                    mpn_mod_sub(res, initial, t, ctx);
                else
                    mpn_mod_add(res, initial, t, ctx);
            }
        }
        return GR_SUCCESS;
    }

    if (n == 2)
    {
        ulong A0, A1, B0, B1;
        ulong p3, p2, p1, p0;
        ulong s4, s3, s2, s1, s0;
        ulong u2, u1;
        ulong v3, v2;

        s4 = s3 = s2 = s1 = s0 = 0;
        u2 = u1 = 0;
        v3 = v2 = 0;

        for (i = 0; i < len; i++)
        {
            A0 = vec1[2 * i + 0];
            A1 = vec1[2 * i + 1];
            B0 = vec2[2 * i + 0];
            B1 = vec2[2 * i + 1];

            umul_ppmm(p2, p1, A1, B0);
            add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), p2, p1);
            umul_ppmm(p1, p0, B0, A0);
            add_sssaaaaaa(u2, u1, s0, u2, u1, s0, UWORD(0), p1, p0);
            umul_ppmm(p2, p1, B1, A0);
            add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), p2, p1);
            umul_ppmm(p3, p2, B1, A1);
            add_sssaaaaaa(s4, v3, v2, s4, v3, v2, UWORD(0), p3, p2);
        }

        /* s3 is small, so this doesn't overflow */
        add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), u2, u1);
        add_sssaaaaaa(s4, s3, s2, s4, s3, s2, UWORD(0), v3, v2);

        s[0] = s0;
        s[1] = s1;
        s[2] = s2;
        s[3] = s3;
        s[4] = s4;
    }
    else
    {
        flint_mpn_mul_n(s, vec1, vec2, n);
        s[2 * n] = 0;
        for (i = 1; i < len; i++)
        {
            flint_mpn_mul_n(t, vec1 + i * n, vec2 + i * n, n);
            s[2 * n] += mpn_add_n(s, s, t, 2 * n);
        }
    }

    sn = 2 * n + (s[2 * n] != 0);

    if (initial == NULL)
    {
        if (!subtract)
        {
            mpn_mod_set_mpn(res, s, sn, ctx);
        }
        else
        {
            mpn_mod_set_mpn(t, s, sn, ctx);
            mpn_mod_neg(res, t, ctx);
        }
    }
    else
    {
        /* todo: consider setting initial as a starting value for the sum */
        mpn_mod_set_mpn(t, s, sn, ctx);

        if (!subtract)
            mpn_mod_add(res, initial, t, ctx);
        else
            mpn_mod_sub(res, initial, t, ctx);
    }

    return GR_SUCCESS;
}

/* todo: optimize for length 1, 2 */
/* todo: optimize for when 2n rather than 2n+1 limbs suffice */
int
_mpn_mod_vec_dot_rev(nn_ptr res, nn_srcptr initial, int subtract, nn_srcptr vec1, nn_srcptr vec2, slong len, gr_ctx_t ctx)
{
    ulong s[2 * MPN_MOD_MAX_LIMBS + 1];
    ulong t[2 * MPN_MOD_MAX_LIMBS];
    slong n = MPN_MOD_CTX_NLIMBS(ctx);
    slong sn;
    slong i;

    if (len <= 2)
    {
        if (len == 0)
        {
            if (initial == NULL)
                flint_mpn_zero(res, n);
            else
                flint_mpn_copyi(res, initial, n);
        }
        else
        {
            if (initial == NULL)
            {
                if (len == 1)
                    mpn_mod_mul(res, vec1, vec2, ctx);
                else
                    mpn_mod_fmma(res, vec1, vec2 + n, vec1 + n, vec2, ctx);

                if (subtract)
                    mpn_mod_neg(res, res, ctx);
            }
            else
            {
                if (len == 1)
                    mpn_mod_mul(t, vec1, vec2, ctx);
                else
                    mpn_mod_fmma(t, vec1, vec2 + n, vec1 + n, vec2, ctx);

                if (subtract)
                    mpn_mod_sub(res, initial, t, ctx);
                else
                    mpn_mod_add(res, initial, t, ctx);
            }
        }
        return GR_SUCCESS;
    }

    if (n == 2)
    {
        ulong A0, A1, B0, B1;
        ulong p3, p2, p1, p0;
        ulong s4, s3, s2, s1, s0;
        ulong u2, u1;
        ulong v3, v2;

        s4 = s3 = s2 = s1 = s0 = 0;
        u2 = u1 = 0;
        v3 = v2 = 0;

        for (i = 0; i < len; i++)
        {
            A0 = vec1[2 * i + 0];
            A1 = vec1[2 * i + 1];
            B0 = vec2[2 * (len - 1 - i) + 0];
            B1 = vec2[2 * (len - 1 - i) + 1];

            umul_ppmm(p2, p1, A1, B0);
            add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), p2, p1);
            umul_ppmm(p1, p0, B0, A0);
            add_sssaaaaaa(u2, u1, s0, u2, u1, s0, UWORD(0), p1, p0);
            umul_ppmm(p2, p1, B1, A0);
            add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), p2, p1);
            umul_ppmm(p3, p2, B1, A1);
            add_sssaaaaaa(s4, v3, v2, s4, v3, v2, UWORD(0), p3, p2);
        }

        /* s3 is small, so this doesn't overflow */
        add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), u2, u1);
        add_sssaaaaaa(s4, s3, s2, s4, s3, s2, UWORD(0), v3, v2);

        s[0] = s0;
        s[1] = s1;
        s[2] = s2;
        s[3] = s3;
        s[4] = s4;
    }
    else
    {
        flint_mpn_mul_n(s, vec1, vec2 + (len - 1) * n, n);
        s[2 * n] = 0;
        for (i = 1; i < len; i++)
        {
            flint_mpn_mul_n(t, vec1 + i * n, vec2 + (len - i - 1) * n, n);
            s[2 * n] += mpn_add_n(s, s, t, 2 * n);
        }
    }

    sn = 2 * n + (s[2 * n] != 0);

    if (initial == NULL)
    {
        if (!subtract)
        {
            mpn_mod_set_mpn(res, s, sn, ctx);
        }
        else
        {
            mpn_mod_set_mpn(t, s, sn, ctx);
            mpn_mod_neg(res, t, ctx);
        }
    }
    else
    {
        /* todo: consider setting initial as a starting value for the sum */
        mpn_mod_set_mpn(t, s, sn, ctx);

        if (!subtract)
            mpn_mod_add(res, initial, t, ctx);
        else
            mpn_mod_sub(res, initial, t, ctx);
    }

    return GR_SUCCESS;
}
