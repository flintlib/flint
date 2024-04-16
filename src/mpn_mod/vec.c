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
_mpn_mod_vec_clear(mp_ptr FLINT_UNUSED(res), slong FLINT_UNUSED(len), gr_ctx_t FLINT_UNUSED(ctx))
{
    return GR_SUCCESS;
}

int
_mpn_mod_vec_zero(mp_ptr res, slong len, gr_ctx_t ctx)
{
    flint_mpn_zero(res, len * MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

int
_mpn_mod_vec_set(mp_ptr res, mp_srcptr x, slong len, gr_ctx_t ctx)
{
    flint_mpn_copyi(res, x, len * MPN_MOD_CTX_NLIMBS(ctx));
    return GR_SUCCESS;
}

void
_mpn_mod_vec_swap(mp_ptr vec1, mp_ptr vec2, slong len, gr_ctx_t ctx)
{
    slong i;
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    for (i = 0; i < len * n; i++)
        FLINT_SWAP(mp_limb_t, vec1[i], vec2[i]);
}

int
_mpn_mod_vec_neg(mp_ptr res, mp_srcptr x, slong len, gr_ctx_t ctx)
{
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    mp_srcptr d = MPN_MOD_CTX_MODULUS(ctx);
    slong i;

    if (n == 2)
    {
        /* Only read to registers once */
        mp_limb_t dd[2];
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
_mpn_mod_vec_add(mp_ptr res, mp_srcptr x, mp_srcptr y, slong len, gr_ctx_t ctx)
{
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    mp_srcptr d = MPN_MOD_CTX_MODULUS(ctx);
    slong i;

    if (n == 2)
    {
        /* Only read to registers once */
        mp_limb_t dd[2];
        dd[0] = d[0];
        dd[1] = d[1];

        for (i = 0; i < len; i++)
            flint_mpn_addmod_2(res + i * n, x + i * n, y + i * n, dd);
    }
    else
        for (i = 0; i < len; i++)
            flint_mpn_addmod_n(res + i * n, x + i * n, y + i * n, d, n);

    return GR_SUCCESS;
}

int
_mpn_mod_vec_sub(mp_ptr res, mp_srcptr x, mp_srcptr y, slong len, gr_ctx_t ctx)
{
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    mp_srcptr d = MPN_MOD_CTX_MODULUS(ctx);
    slong i;

    if (n == 2)
    {
        /* Only read to registers once */
        mp_limb_t dd[2];
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
_mpn_mod_vec_mul(mp_ptr res, mp_srcptr x, mp_srcptr y, slong len, gr_ctx_t ctx)
{
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    slong i;

    if (n == 2)
    {
        mp_bitcnt_t norm = MPN_MOD_CTX_NORM(ctx);
        mp_srcptr dnormed = MPN_MOD_CTX_MODULUS_NORMED(ctx);
        mp_srcptr dinv = MPN_MOD_CTX_MODULUS_PREINV(ctx);

        for (i = 0; i < len; i++)
            flint_mpn_mulmod_preinvn_2(res + i * n, x + i * n, y + i * n, dnormed, dinv, norm);
    }
    else
        for (i = 0; i < len; i++)
            mpn_mod_mul(res + i * n, x + i * n, y + i * n, ctx);

    return GR_SUCCESS;
}

/* todo: worth it to check for special cases (0, 1)? */
/* todo: shoup multiplication */
int
_mpn_mod_vec_mul_scalar(mp_ptr res, mp_srcptr x, slong len, mp_srcptr y, gr_ctx_t ctx)
{
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    slong i;

    if (n == 2)
    {
        mp_bitcnt_t norm = MPN_MOD_CTX_NORM(ctx);
        mp_srcptr dnormed = MPN_MOD_CTX_MODULUS_NORMED(ctx);
        mp_srcptr dinv = MPN_MOD_CTX_MODULUS_PREINV(ctx);

        for (i = 0; i < len; i++)
            flint_mpn_mulmod_preinvn_2(res + i * n, x + i * n, y, dnormed, dinv, norm);
    }
    else
        for (i = 0; i < len; i++)
            mpn_mod_mul(res + i * n, x + i * n, y, ctx);

    return GR_SUCCESS;
}

int
_mpn_mod_scalar_mul_vec(mp_ptr res, mp_srcptr y, mp_srcptr x, slong len, gr_ctx_t ctx)
{
    return _mpn_mod_vec_mul_scalar(res, x, len, y, ctx);
}

/* todo: worth it to check for special cases (0, 1)? */
/* todo: shoup multiplication */
int
_mpn_mod_vec_addmul_scalar(mp_ptr res, mp_srcptr x, slong len, mp_srcptr y, gr_ctx_t ctx)
{
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    slong i;

    if (n == 2)
    {
        mp_limb_t t[2];
        mp_bitcnt_t norm = MPN_MOD_CTX_NORM(ctx);
        mp_srcptr dnormed = MPN_MOD_CTX_MODULUS_NORMED(ctx);
        mp_srcptr dinv = MPN_MOD_CTX_MODULUS_PREINV(ctx);
        mp_srcptr d = MPN_MOD_CTX_MODULUS(ctx);

        for (i = 0; i < len; i++)
        {
            flint_mpn_mulmod_preinvn_2(t, x + i * n, y, dnormed, dinv, norm);
            flint_mpn_addmod_2(res + i * n, res + i * n, t, d);
        }
    }
    else
    {
        mp_limb_t t[MPN_MOD_MAX_LIMBS];

        for (i = 0; i < len; i++)
        {
            mpn_mod_mul(t, x + i * n, y, ctx);
            mpn_mod_add(res + i * n, res + i * n, t, ctx);
        }
    }

    return GR_SUCCESS;
}

/* todo: optimize for length 1, 2 */
/* todo: optimize for when 2n rather than 2n+1 limbs suffice */
int
_mpn_mod_vec_dot(mp_ptr res, mp_srcptr initial, int subtract, mp_srcptr vec1, mp_srcptr vec2, slong len, gr_ctx_t ctx)
{
    mp_limb_t s[2 * MPN_MOD_MAX_LIMBS + 1];
    mp_limb_t t[2 * MPN_MOD_MAX_LIMBS];
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    mp_size_t sn;
    slong i;

    if (len <= 0)
    {
        if (initial == NULL)
            flint_mpn_zero(res, n);
        else
            flint_mpn_copyi(res, initial, n);
        return GR_SUCCESS;
    }

    if (n == 2)
    {
        mp_limb_t A0, A1, B0, B1;
        mp_limb_t p3, p2, p1, p0;
        mp_limb_t s4, s3, s2, s1, s0;
        mp_limb_t u2, u1;
        mp_limb_t v3, v2;

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
_mpn_mod_vec_dot_rev(mp_ptr res, mp_srcptr initial, int subtract, mp_srcptr vec1, mp_srcptr vec2, slong len, gr_ctx_t ctx)
{
    mp_limb_t s[2 * MPN_MOD_MAX_LIMBS + 1];
    mp_limb_t t[2 * MPN_MOD_MAX_LIMBS];
    mp_size_t n = MPN_MOD_CTX_NLIMBS(ctx);
    mp_size_t sn;
    slong i;

    if (len <= 0)
    {
        if (initial == NULL)
            flint_mpn_zero(res, n);
        else
            flint_mpn_copyi(res, initial, n);
        return GR_SUCCESS;
    }

    if (n == 2)
    {
        mp_limb_t A0, A1, B0, B1;
        mp_limb_t p3, p2, p1, p0;
        mp_limb_t s4, s3, s2, s1, s0;
        mp_limb_t u2, u1;
        mp_limb_t v3, v2;

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
