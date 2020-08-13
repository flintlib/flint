/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"
#include "nmod_vec.h"


void n_fq_print_pretty(
    const mp_limb_t * a,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    int first;

    first = 1;
    for (i = d - 1; i >= 0; i--)
    {
        if (a[i] == 0)
            continue;

        if (!first)
            flint_printf("+");

        first = 0;
        flint_printf("%wu", a[i]);

        if (i > 0)
        {
            flint_printf("*%s", ctx->var);
            if (i > 1)
                flint_printf("^%wd", i);
        }
    }

    if (first)
        flint_printf("0");
}


void n_fq_get_fq_nmod(
    fq_nmod_t a,
    const mp_limb_t * b,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    slong d = fq_nmod_ctx_degree(ctx);

    nmod_poly_fit_length(a, d);

    for (i = 0; i < d; i++)
        a->coeffs[i] = b[i];

    a->length = d;
    _nmod_poly_normalise(a);
}

void n_fq_set_fq_nmod(
    mp_limb_t * a,
    const fq_nmod_t b,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    slong d = fq_nmod_ctx_degree(ctx);

    FLINT_ASSERT(b->length <= d);

    for (i = 0; i < d; i++)
        a[i] = i < b->length ? b->coeffs[i] : 0;
}

void n_fq_add_si(
    mp_limb_t * a,
    const mp_limb_t * b,
    slong c,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);

    if (a != b)
        _nmod_vec_set(a, b, d);

    if (c < 0)
    {
        ulong cc = -c;
        if (cc >= ctx->mod.n)
            NMOD_RED(cc, cc, ctx->mod);
        a[0] = nmod_sub(a[0], cc, ctx->mod);
    }
    else
    {
        ulong cc = c;
        if (cc >= ctx->mod.n)
            NMOD_RED(cc, cc, ctx->mod);
        a[0] = nmod_add(a[0], cc, ctx->mod);
    }
}

void _n_fq_reduce(
    mp_limb_t * a,
    mp_limb_t * b,
    slong blen,
    const fq_nmod_ctx_t ctx,
    mp_limb_t * t)  /* length 2d */
{
    slong i, j, k, deg = ctx->modulus->length - 1;
    slong d = ctx->j[ctx->len - 1];

    FLINT_ASSERT(a != b);

    FLINT_ASSERT(0 <= blen && blen <= 2*d - 1);
    FLINT_ASSERT(blen == 0 || b[blen - 1] != 0);

    if (blen <= deg)
    {
        for (i = 0; i < blen; i++)
            a[i] = b[i];
        for (i = blen; i < deg; i++)
            a[i] = 0;
    }
    else if (ctx->sparse_modulus)
    {
        nmod_t mod = ctx->mod;

        for (k = ctx->len - 2; k >= 0; k--)
            t[k] = mod.n - ctx->a[k];

        for (i = blen - 1; i >= d; i--)
        {
            for (k = ctx->len - 2; k >= 0; k--)
                NMOD_ADDMUL(b[ctx->j[k] + i - d], b[i], t[k], mod);

            b[i] = 0;
        }

        for (i = 0; i < deg; i++)
            a[i] = b[i];
    }
    else
    {
/*
        _nmod_poly_divrem_newton_n_preinv(t, a, b, blen, 
                                          ctx->modulus->coeffs, ctx->modulus->length,
                                          ctx->inv->coeffs, ctx->inv->length,
                                          ctx->mod);
*/
        mp_limb_t * Q = t;
        mp_limb_t * R = a;
        const mp_limb_t * A = b;
        slong lenA = blen;
        const mp_limb_t * B = ctx->modulus->coeffs;
        slong lenB = deg + 1;
        const mp_limb_t * Binv = ctx->inv->coeffs;
        slong lenBinv = ctx->inv->length;

        const slong lenQ = lenA - lenB + 1;

        FLINT_ASSERT(lenBinv > 0);
        FLINT_ASSERT(lenQ > 0);

        if (lenQ <= 2)
        {
            if (lenQ == 2)
                _nmod_poly_divrem_q1(Q, R, A, lenA, B, lenB, ctx->mod);
            else
                _nmod_poly_divrem_q0(Q, R, A, B, lenB, ctx->mod);
            return;
        }

        if (deg < 20)
        {
            for (i = 0; i < lenQ; i++)
            {
                mp_limb_t t2 = 0, t1 = 0, t0 = 0;
                j = FLINT_MAX(0, i - lenBinv + 1);
                umul_ppmm(t1, t0, A[lenA - 1 - j], Binv[i - j]);
                for (j++; j <= i; j++)
                    MAC(t2, t1, t0, A[lenA - 1 - j], Binv[i - j]);
                NMOD_RED3(Q[lenQ - 1 - i], t2, t1, t0, ctx->mod);
            }

            for (i = 0; i < deg; i++)
            {
                mp_limb_t t2 = 0, t1 = 0, t0 = 0;
                for (j = FLINT_MAX(0, i - lenQ + 1); j <= i; j++)
                    MAC(t2, t1, t0, B[j], Q[i - j]);
                NMOD_RED3(t0, t2, t1, t0, ctx->mod);
                R[i] = nmod_sub(A[i], t0, ctx->mod);
            }
        }
        else
        {
            mp_ptr Arev = t + d;

            _nmod_poly_reverse(Arev, A + (lenA - lenQ), lenQ, lenQ);

            _nmod_poly_mullow(Q, Arev, lenQ, Binv, FLINT_MIN(lenQ, lenBinv), lenQ, ctx->mod);

            _nmod_poly_reverse(Q, Q, lenQ, lenQ);

            FLINT_ASSERT(lenB > 1);
            FLINT_ASSERT(lenQ < lenB - 1);

            _nmod_poly_mullow(R, B, lenB - 1, Q, lenQ, lenB - 1, ctx->mod);

            _nmod_vec_sub(R, A, R, lenB - 1, ctx->mod);
        }
    }
}

void _n_fq_madd2(
    mp_limb_t * a,          /* length 2d */
    const mp_limb_t * b,    /* length d */
    const mp_limb_t * c,    /* length d */
    const fq_nmod_ctx_t ctx,
    mp_limb_t * t)          /* length 2d */
{
    slong d = ctx->modulus->length - 1;
    FLINT_ASSERT(d > 0);
    if (d < 30)
    {
        slong i, j;
        for (i = 0; i + 1 < d; i++)
        {
            mp_limb_t t2 = 0, t1 = 0, t0 = 0;
            mp_limb_t s2 = 0, s1 = 0, s0 = 0;
            umul_ppmm(t1, t0, b[i], c[0]);
            umul_ppmm(s1, s0, b[d - 1], c[d - 1 - i]);

            add_ssaaaa(t1, t0, t1, t0, 0, a[i]);
            add_ssaaaa(s1, s0, s1, s0, 0, a[2*d - 2 - i]);

            for (j = 1; j <= i; j++)
            {
                MAC(t2, t1, t0, b[i - j], c[0 + j]);
                MAC(s2, s1, s0, b[d - 1 - j], c[d - 1 - i + j]);
            }
            NMOD_RED3(a[i], t2, t1, t0, ctx->mod);
            NMOD_RED3(a[2*d - 2 - i], s2, s1, s0, ctx->mod);
        }

        {
            mp_limb_t t2 = 0, t1 = 0, t0 = 0;
            umul_ppmm(t1, t0, b[d - 1], c[0]);
            add_ssaaaa(t1, t0, t1, t0, 0, a[d - 1]);
            for (j = 1; j < d; j++)
            {
                MAC(t2, t1, t0, b[d - 1 - j], c[0 + j]);
            }
            NMOD_RED3(a[d - 1], t2, t1, t0, ctx->mod);
        }
    }
    else
    {
        _nmod_poly_mul(t, b, d, c, d, ctx->mod);
        _nmod_vec_add(a, a, t, 2*d - 1, ctx->mod);
    }
}

void _n_fq_mul2(
    mp_limb_t * a,          /* length 2d */
    const mp_limb_t * b,    /* length d */
    const mp_limb_t * c,    /* length d */
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    FLINT_ASSERT(d > 0);
    if (d < 30)
    {
        slong i, j;
        for (i = 0; i + 1 < d; i++)
        {
            mp_limb_t t2 = 0, t1 = 0, t0 = 0;
            mp_limb_t s2 = 0, s1 = 0, s0 = 0;
            umul_ppmm(t1, t0, b[i], c[0]);
            umul_ppmm(s1, s0, b[d - 1], c[d - 1 - i]);
            for (j = 1; j <= i; j++)
            {
                MAC(t2, t1, t0, b[i - j], c[0 + j]);
                MAC(s2, s1, s0, b[d - 1 - j], c[d - 1 - i + j]);
            }
            NMOD_RED3(a[i], t2, t1, t0, ctx->mod);
            NMOD_RED3(a[2*d - 2 - i], s2, s1, s0, ctx->mod);
        }

        {
            mp_limb_t t2 = 0, t1 = 0, t0 = 0;
            umul_ppmm(t1, t0, b[d - 1], c[0]);
            for (j = 1; j < d; j++)
            {
                MAC(t2, t1, t0, b[d - 1 - j], c[0 + j]);
            }
            NMOD_RED3(a[d - 1], t2, t1, t0, ctx->mod);
        }
    }
    else
    {
        _nmod_poly_mul(a, b, d, c, d, ctx->mod);
    }
}

void _n_fq_inv(
    mp_limb_t * a,
    const mp_limb_t * b,
    const fq_nmod_ctx_t ctx,
    mp_limb_t * t) /* length d */
{
    slong d = ctx->modulus->length - 1;
    slong blen = d;
    FLINT_ASSERT(d > 0);

    while (blen > 0 && b[blen - 1] == 0)
        blen--;

    if (blen < 1)
    {
        flint_throw(FLINT_ERROR, "impossible inverse in _fq_nmod_inv");
    }
    else if (blen == 1)
    {
        a[0] = n_invmod(b[0], ctx->mod.n);
        _nmod_vec_zero(a + 1, d - 1);
    }
    else
    {
        if (1 != _nmod_poly_gcdinv(t, a, b, blen, ctx->modulus->coeffs, d + 1, ctx->mod))
        {
            flint_throw(FLINT_ERROR, "impossible inverse in _fq_nmod_inv");
        }

        if (t[0] != 1)
        {
            _nmod_vec_scalar_mul_nmod(a, a, d, n_invmod(t[0], ctx->mod.n), ctx->mod);
        }
    }
}

void n_fq_mul(
    mp_limb_t * a,
    const mp_limb_t * b,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx)
{
    fq_nmod_t A, B, C;
    fq_nmod_init(A, ctx);
    fq_nmod_init(B, ctx);
    fq_nmod_init(C, ctx);
    n_fq_get_fq_nmod(B, b, ctx);
    n_fq_get_fq_nmod(C, c, ctx);
    fq_nmod_mul(A, B, C, ctx);
    n_fq_set_fq_nmod(a, A, ctx);
    fq_nmod_clear(A, ctx);
    fq_nmod_clear(B, ctx);
    fq_nmod_clear(C, ctx);
}

void n_fq_inv(
    mp_limb_t * a,
    const mp_limb_t * b,
    const fq_nmod_ctx_t ctx)
{
    fq_nmod_t A, B;
    fq_nmod_init(A, ctx);
    fq_nmod_init(B, ctx);
    n_fq_get_fq_nmod(B, b, ctx);
    fq_nmod_inv(A, B, ctx);
    n_fq_set_fq_nmod(a, A, ctx);
    fq_nmod_clear(A, ctx);
    fq_nmod_clear(B, ctx);
}

void _n_fq_pow_ui(
    mp_limb_t * a,
    const mp_limb_t * b,
    ulong e,
    const fq_nmod_ctx_t ctx)
{
    fq_nmod_t A, B;
    fq_nmod_init(A, ctx);
    fq_nmod_init(B, ctx);
    n_fq_get_fq_nmod(B, b, ctx);
    fq_nmod_pow_ui(A, B, e, ctx);
    n_fq_set_fq_nmod(a, A, ctx);
    fq_nmod_clear(A, ctx);
    fq_nmod_clear(B, ctx);
}

void n_fq_pow_ui(
    mp_limb_t * a,
    const mp_limb_t * b,
    ulong e,
    const fq_nmod_ctx_t ctx)
{
    fq_nmod_t A, B;
    fq_nmod_init(A, ctx);
    fq_nmod_init(B, ctx);
    n_fq_get_fq_nmod(B, b, ctx);
    fq_nmod_pow_ui(A, B, e, ctx);
    n_fq_set_fq_nmod(a, A, ctx);
    fq_nmod_clear(A, ctx);
    fq_nmod_clear(B, ctx);
}

int n_fq_is_canonical(
    const mp_limb_t * a,
    const fq_nmod_ctx_t ctx)
{
    slong i, d = fq_nmod_ctx_degree(ctx);
    for (i = 0; i < d; i++)
    {
        if (a[i] >= ctx->mod.n)
            return 0;
    }
    return 1;
}
