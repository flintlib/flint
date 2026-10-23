/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "mpn_extras.h"
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_poly.h"
#include "fq_nmod.h"
#include "fq_zech.h"
#include "fq_zech_poly.h"

/* Writes the d coordinates of x (with respect to the power basis) to res,
   given the magic multiplier pinv for division by p (see below). */
static void
_fq_zech_get_coeffs(ulong * res, ulong x, slong d, ulong p, ulong pinv,
    const fq_zech_ctx_struct * ctx)
{
    ulong q, r, hi, lo;
    slong j;

    if (x == ctx->qm1)
    {
        for (j = 0; j < d; j++)
            res[j] = 0;
        return;
    }

    q = ctx->eval_table[x];

    if (p == 2)
    {
        for (j = 0; j < d; j++)
        {
            res[j] = q & 1;
            q >>= 1;
        }
    }
    else if (pinv != 0)
    {
        for (j = 0; j < d; j++)
        {
            umul_ppmm(hi, lo, q, pinv);
            (void) lo;
            res[j] = q - hi * p;
            q = hi;
        }
    }
    else
    {
        for (j = 0; j < d; j++)
        {
            r = n_divrem2_precomp(&q, q, p, ctx->ppre);
            res[j] = r;
        }
    }
}

/* Zech logarithm of sum_j c[j] * a^j where 0 <= c[j] < p. */
static ulong
_fq_zech_log_from_coeffs(const ulong * c, slong d, const fq_zech_ctx_struct * ctx)
{
    ulong qm1 = ctx->qm1;
    const ulong * table = ctx->zech_log_table;
    const ulong * pft = ctx->prime_field_table;
    ulong s0, s1, t, u;
    slong j;

    /* Two independent chains of table lookups to hide latency. */
    s0 = s1 = qm1;

    for (j = 0; j + 1 < d; j += 2)
    {
        if (c[j] != 0)
        {
            t = pft[c[j]] + j;
            /* no wraparound: the log of a prime field element is a
               multiple k (q - 1) / (p - 1) with k <= p - 2, and j < d */
            FLINT_ASSERT(t < qm1);

            if (s0 == qm1)
                s0 = t;
            else
            {
                u = table[n_submod(s0, t, qm1)];
                s0 = (u == qm1) ? qm1 : n_addmod(u, t, qm1);
            }
        }

        if (c[j + 1] != 0)
        {
            t = pft[c[j + 1]] + j + 1;
            FLINT_ASSERT(t < qm1);

            if (s1 == qm1)
                s1 = t;
            else
            {
                u = table[n_submod(s1, t, qm1)];
                s1 = (u == qm1) ? qm1 : n_addmod(u, t, qm1);
            }
        }
    }

    if (j < d && c[j] != 0)
    {
        t = pft[c[j]] + j;
        FLINT_ASSERT(t < qm1);

        if (s0 == qm1)
            s0 = t;
        else
        {
            u = table[n_submod(s0, t, qm1)];
            s0 = (u == qm1) ? qm1 : n_addmod(u, t, qm1);
        }
    }

    if (s0 == qm1)
        return s1;
    if (s1 == qm1)
        return s0;

    u = table[n_submod(s0, s1, qm1)];
    return (u == qm1) ? qm1 : n_addmod(u, s1, qm1);
}

/* Reduces R (of length 2d - 1, entries in [0, p)) modulo the defining
   polynomial, lazily accumulating when (d - 1)(p - 1)^2 + p fits in a word. */
static void
_fq_zech_reduce_coeffs(ulong * R, slong d, int lazy,
    const fq_nmod_ctx_struct * fctx)
{
    const ulong * a = fctx->a;
    const slong * jj = fctx->j;
    const slong len = fctx->len;
    const ulong p = fctx->mod.n;
    slong i, k;
    ulong c;

    if (d == 1)
        return;

    if (p == 2)
    {
        for (i = 2 * d - 2; i >= d; i--)
        {
            c = R[i];
            if (c != 0)
                for (k = 0; k < len - 1; k++)
                    R[jj[k] + i - d] ^= 1;
        }
    }
    else if (lazy)
    {
        for (i = 2 * d - 2; i >= d; i--)
        {
            NMOD_RED(c, R[i], fctx->mod);
            if (c != 0)
                for (k = 0; k < len - 1; k++)
                    R[jj[k] + i - d] += c * (p - a[k]);
        }

        for (i = 0; i < d; i++)
            NMOD_RED(R[i], R[i], fctx->mod);
    }
    else
    {
        _fq_nmod_reduce(R, 2 * d - 1, fctx);
    }
}

/* Computes coefficients [nlo, nhi) of the product by packing the
   coordinates of the coefficients into polynomials over Z/pZ with stride
   2d - 1 (Kronecker substitution), multiplying over Z/pZ, and reducing
   each block of 2d - 1 coordinates modulo the defining polynomial. */
void
_fq_zech_poly_mulmid_univariate(fq_zech_struct * rop,
    const fq_zech_struct * op1, slong len1,
    const fq_zech_struct * op2, slong len2,
    slong nlo, slong nhi, const fq_zech_ctx_t ctx)
{
    const fq_nmod_ctx_struct * fctx = ctx->fq_nmod_ctx;
    const slong d = fq_zech_ctx_degree(ctx);
    const slong pd = 2 * d - 1;
    const ulong qm1 = ctx->qm1;
    const ulong p = ctx->p;
    ulong pinv, hi, lo, v;
    slong i, j, m, clen1, clen2, cm, trunc;
    int squaring, lazy;
    nn_ptr cop1, cop2, crop, rev;

    squaring = (op1 == op2 && len1 == len2);

    len1 = FLINT_MIN(len1, nhi);
    len2 = FLINT_MIN(len2, nhi);

    while (len1 > 0 && op1[len1 - 1].value == qm1)
        len1--;
    while (len2 > 0 && op2[len2 - 1].value == qm1)
        len2--;

    m = (len1 == 0 || len2 == 0) ? 0 : FLINT_MIN(nhi, len1 + len2 - 1);

    if (m <= nlo)
    {
        for (i = 0; i < nhi - nlo; i++)
            rop[i].value = qm1;
        return;
    }

    /* Low coefficients of either input that only contribute to
       coefficients below nlo can be discarded. */
    if (nlo != 0 && !squaring)
    {
        trunc = len1 - (len1 + len2 - 1 - nlo);
        if (trunc > 0)
        {
            op1 += trunc;
            len1 -= trunc;
            nlo -= trunc;
            nhi -= trunc;
            m -= trunc;
        }

        trunc = len2 - (len1 + len2 - 1 - nlo);
        if (trunc > 0)
        {
            op2 += trunc;
            len2 -= trunc;
            nlo -= trunc;
            nhi -= trunc;
            m -= trunc;
        }
    }

    /* floor(x / p) = floor(x * pinv / 2^FLINT_BITS) for x, p < 2^32 */
#if FLINT_BITS == 64
    pinv = (qm1 < (UWORD(1) << 32)) ? (UWORD_MAX / p + 1) : 0;
#else
    pinv = 0;
#endif

    clen1 = pd * (len1 - 1) + d;
    clen2 = pd * (len2 - 1) + d;
    /* clen1 + clen2 - 1 = pd * (len1 + len2 - 1) >= pd * m */
    cm = pd * m;
    FLINT_ASSERT(cm <= clen1 + clen2 - 1);

    cop1 = flint_malloc(sizeof(ulong) *
        (clen1 + (squaring ? 0 : clen2) + pd * (m - nlo)));
    cop2 = squaring ? cop1 : cop1 + clen1;
    crop = cop2 + (squaring ? clen1 : clen2);

    for (i = 0; i < len1; i++)
    {
        _fq_zech_get_coeffs(cop1 + pd * i, op1[i].value, d, p, pinv, ctx);
        if (i < len1 - 1)
            flint_mpn_zero(cop1 + pd * i + d, d - 1);
    }

    if (!squaring)
    {
        for (i = 0; i < len2; i++)
        {
            _fq_zech_get_coeffs(cop2 + pd * i, op2[i].value, d, p, pinv, ctx);
            if (i < len2 - 1)
                flint_mpn_zero(cop2 + pd * i + d, d - 1);
        }
    }

    if (nlo == 0)
    {
        if (clen1 >= clen2)
            _nmod_poly_mullow(crop, cop1, clen1, cop2, clen2, cm, fctx->mod);
        else
            _nmod_poly_mullow(crop, cop2, clen2, cop1, clen1, cm, fctx->mod);
    }
    else
    {
        _nmod_poly_mulmid(crop, cop1, clen1, cop2, clen2, pd * nlo, cm,
                          fctx->mod);
    }

    /* (d - 1)(p - 1)^2 + p must fit in a word for lazy reduction */
    umul_ppmm(hi, lo, p - 1, p - 1);
    lazy = (hi == 0);
    if (lazy)
    {
        umul_ppmm(hi, lo, lo, d - 1);
        lazy = (hi == 0) && (lo + p >= lo);
    }

    /* When there are many outputs, it pays off to build the table of
       discrete logarithms (the inverse of the evaluation table) so that
       the logarithm of each output coefficient costs a single lookup. */
    if ((ulong) (m - nlo) * d >= qm1)
    {
        rev = flint_malloc(sizeof(ulong) * (qm1 + 1));

        for (v = 0; v <= qm1; v++)
            rev[ctx->eval_table[v]] = v;

        for (i = 0; i < m - nlo; i++)
        {
            nn_ptr R = crop + pd * i;

            _fq_zech_reduce_coeffs(R, d, lazy, fctx);

            v = R[d - 1];
            if (p == 2)
                for (j = d - 2; j >= 0; j--)
                    v = 2 * v + R[j];
            else
                for (j = d - 2; j >= 0; j--)
                    v = v * p + R[j];

            rop[i].value = rev[v];
        }

        flint_free(rev);
    }
    else
    {
        for (i = 0; i < m - nlo; i++)
        {
            _fq_zech_reduce_coeffs(crop + pd * i, d, lazy, fctx);
            rop[i].value = _fq_zech_log_from_coeffs(crop + pd * i, d, ctx);
        }
    }

    for ( ; i < nhi - nlo; i++)
        rop[i].value = qm1;

    flint_free(cop1);
}

/* Number of terms a_i b_j with i + j < n in the product of polynomials
   of length len1 and len2. */
static double
_mul_work(slong len1, slong len2, slong n)
{
    double a, b, t;

    a = FLINT_MIN(FLINT_MIN(len1, len2), n);
    b = FLINT_MIN(FLINT_MAX(len1, len2), n);

    if (n <= a)
        return (double) n * (n + 1) / 2;
    else if (n <= b)
        return a * (a + 1) / 2 + (n - a) * a;

    t = a + b - 1 - n;
    t = FLINT_MAX(t, 0);
    return a * b - t * (t + 1) / 2;
}

/* Crossover (average number of terms per output coefficient) between
   classical and univariate multiplication. The classical algorithm costs
   one Zech logarithm table lookup per term; the univariate algorithm costs
   O(d) per output coefficient (packing, multiplication over Z/pZ,
   reduction, and O(d) table lookups for the conversion back). Lookups get
   more expensive when the tables do not fit in cache (the tables take
   2q words), making univariate multiplication competitive earlier. */
slong
_fq_zech_poly_mul_univariate_threshold(const fq_zech_ctx_t ctx)
{
    slong d = fq_zech_ctx_degree(ctx);

    if (ctx->qm1 < (UWORD(1) << 17))
        return 9 * d - 4;
    else if (ctx->qm1 < (UWORD(1) << 18))
        return FLINT_MAX((9 * d - 4) / 2, 2 * d + 2);
    else
        return 2 * d + 2;
}

int
_fq_zech_poly_mulmid_want_univariate(slong len1, slong len2,
    slong nlo, slong nhi, const fq_zech_ctx_t ctx)
{
    slong T, trunc, nhi2;
    double W;

    len1 = FLINT_MIN(len1, nhi);
    len2 = FLINT_MIN(len2, nhi);
    nhi = FLINT_MIN(nhi, len1 + len2 - 1);

    if (nlo >= nhi)
        return 0;

    T = _fq_zech_poly_mul_univariate_threshold(ctx);

    /* There are at most min(len1, len2) terms per output coefficient. */
    if (FLINT_MIN(len1, len2) < T)
        return 0;

    /* Scaled by 5. Truncated products have a slightly higher crossover
       since the classical algorithm saves more work. */
    if (nlo == 0 && nhi == len1 + len2 - 1)
        T = 5 * T;
    else
        T = 6 * T;

    W = _mul_work(len1, len2, nhi) - _mul_work(len1, len2, nlo);

    /* Length of the Kronecker product after discarding low input terms
       which do not contribute. */
    nhi2 = nhi;
    trunc = nlo - (len2 - 1);
    if (trunc > 0)
        nhi2 -= trunc;
    trunc = nlo - (len1 - 1);
    if (trunc > 0)
        nhi2 -= trunc;

    return 5 * W >= (double) T * nhi2;
}

int
_fq_zech_poly_sqr_want_univariate(slong len, const fq_zech_ctx_t ctx)
{
    slong T = _fq_zech_poly_mul_univariate_threshold(ctx);

    if (len < 2 * T)
        return 0;

    /* len^2 >= (11/10) T (2 len - 1) */
    return 10 * (double) len * len >= 11 * (double) T * (2 * len - 1);
}

void
_fq_zech_poly_mullow_univariate(fq_zech_struct * rop,
    const fq_zech_struct * op1, slong len1,
    const fq_zech_struct * op2, slong len2,
    slong n, const fq_zech_ctx_t ctx)
{
    _fq_zech_poly_mulmid_univariate(rop, op1, len1, op2, len2, 0, n, ctx);
}

void
_fq_zech_poly_mul_univariate(fq_zech_struct * rop,
    const fq_zech_struct * op1, slong len1,
    const fq_zech_struct * op2, slong len2, const fq_zech_ctx_t ctx)
{
    _fq_zech_poly_mulmid_univariate(rop, op1, len1, op2, len2,
                                    0, len1 + len2 - 1, ctx);
}

void
fq_zech_poly_mullow_univariate(fq_zech_poly_t rop,
    const fq_zech_poly_t op1, const fq_zech_poly_t op2,
    slong n, const fq_zech_ctx_t ctx)
{
    const slong len1 = op1->length;
    const slong len2 = op2->length;

    if (len1 == 0 || len2 == 0 || n == 0)
    {
        fq_zech_poly_zero(rop, ctx);
        return;
    }

    n = FLINT_MIN(n, len1 + len2 - 1);

    fq_zech_poly_fit_length(rop, n, ctx);
    _fq_zech_poly_mulmid_univariate(rop->coeffs, op1->coeffs, len1,
                                    op2->coeffs, len2, 0, n, ctx);
    _fq_zech_poly_set_length(rop, n, ctx);
    _fq_zech_poly_normalise(rop, ctx);
}

void
fq_zech_poly_mul_univariate(fq_zech_poly_t rop,
    const fq_zech_poly_t op1, const fq_zech_poly_t op2,
    const fq_zech_ctx_t ctx)
{
    const slong len1 = op1->length;
    const slong len2 = op2->length;

    if (len1 == 0 || len2 == 0)
    {
        fq_zech_poly_zero(rop, ctx);
        return;
    }

    fq_zech_poly_fit_length(rop, len1 + len2 - 1, ctx);
    _fq_zech_poly_mulmid_univariate(rop->coeffs, op1->coeffs, len1,
                                    op2->coeffs, len2, 0, len1 + len2 - 1, ctx);
    _fq_zech_poly_set_length(rop, len1 + len2 - 1, ctx);
    _fq_zech_poly_normalise(rop, ctx);
}
