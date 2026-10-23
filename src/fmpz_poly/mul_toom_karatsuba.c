/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

/* The largest product length worth handling by evaluation/interpolation
   (the Toom scalar code) at a given coefficient size: the O(N^2)
   additions and small multiplications/divisions of the evaluation and
   interpolation are only negligible against the N multiplications for
   large coefficients (measured: below about 2500 bits, splitting with
   Karatsuba down to 3 x 3 is best; N = 15 wins from about 4000 bits). */
static slong
_toom_scalar_nmax(slong bits)
{
    if (bits < 2500)
        return 5;
    else if (bits < 4000)
        return 11;
    else if (bits < 6000)
        return 15;
    else
        return FMPZ_POLY_TOOM_SCALAR_N_MAX;
}

static slong
_toom_scalar_nmax2(const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
{
    slong bits1, bits2;

    bits1 = _fmpz_vec_max_bits(poly1, len1);
    bits1 = FLINT_ABS(bits1);

    if (poly1 == poly2 && len1 == len2)
        bits2 = bits1;
    else
    {
        bits2 = _fmpz_vec_max_bits(poly2, len2);
        bits2 = FLINT_ABS(bits2);
    }

    return _toom_scalar_nmax(FLINT_MAX(bits1, bits2));
}

/*
    Karatsuba multiplication with _fmpz_poly_mul_toom_scalar and
    _fmpz_poly_mulmid_toom_scalar as the base cases, following the
    structure of the mpn-level _flint_mpn_poly_mulmid_karatsuba: full products
    by Karatsuba with chunking for unbalanced lengths, low products by
    Mulders' algorithm, high products by reversal, and middle windows by
    blocks of the Hanrot-Quercia-Zimmermann Karatsuba middle product.

    The base cases are entered by simply attempting the toom_scalar
    calls, which fail (leaving the output untouched) exactly when the
    shape is too large for the tables; every recursion strictly reduces
    the relevant size, so this terminates without duplicating the size
    limit here.

    Intended for coefficients large enough that fmpz arithmetic overhead
    is negligible (roughly 2000 bits and up); the caller is responsible
    for choosing this algorithm only in that regime.
*/

/* Full product, any lengths >= 1. */
static void
_kar_mul(fmpz * res, const fmpz * f, slong flen, const fmpz * g, slong glen)
{
    slong m, f1len, g1len, tlen, ulen, vlen;
    fmpz * t, * u, * v;
    int squaring;

    if (flen < glen)
    {
        FLINT_SWAP(const fmpz *, f, g);
        FLINT_SWAP(slong, flen, glen);
    }

    if (_fmpz_poly_mul_toom_scalar_nmax(res, f, flen, g, glen, _toom_scalar_nmax2(f, flen, g, glen)))
        return;

    squaring = (f == g && flen == glen);

    /* Unbalanced lengths: chunk the longer operand into pieces of the
       length of the shorter one. */
    if (flen > glen + glen / 2)
    {
        slong i, m2, slen = glen;

        t = _fmpz_vec_init(2 * slen - 1);

        _kar_mul(res, f, slen, g, slen);

        for (i = slen; i < flen; i += slen)
        {
            m2 = FLINT_MIN(slen, flen - i);
            _kar_mul(t, f + i, m2, g, slen);
            /* the first slen - 1 coefficients overlap the previous chunk */
            _fmpz_vec_add(res + i, res + i, t, slen - 1);
            _fmpz_vec_set(res + i + slen - 1, t + slen - 1, m2);
        }

        _fmpz_vec_clear(t, 2 * slen - 1);
        return;
    }

    /* split at X = x^m */
    /* res = f0 g0 + ((f0 + f1) (g0 + g1) - f0 g0 - f1 g1) X + f1 g1 X^2 */
    m = (glen + 1) / 2;
    f1len = flen - m;
    g1len = glen - m;

    /* low part: res[0, ..., 2m-2] = f0 g0 */
    _kar_mul(res, f, m, g, m);
    fmpz_zero(res + 2 * m - 1);
    /* high part: res[2m, ..., flen+glen-2] = f1 g1 */
    _kar_mul(res + 2 * m, f + m, f1len, g + m, g1len);

    tlen = FLINT_MAX(m, f1len);
    ulen = FLINT_MAX(m, g1len);
    vlen = tlen + ulen - 1;

    t = _fmpz_vec_init(tlen + (squaring ? 0 : ulen) + vlen);
    v = t + tlen + (squaring ? 0 : ulen);

    _fmpz_poly_add(t, f, m, f + m, f1len);

    if (squaring)
    {
        _kar_mul(v, t, tlen, t, tlen);
    }
    else
    {
        u = t + tlen;
        _fmpz_poly_add(u, g, m, g + m, g1len);
        _kar_mul(v, t, tlen, u, ulen);
    }

    /* v -= f0 g0, v -= f1 g1, res += v X^m */
    _fmpz_vec_sub(v, v, res, 2 * m - 1);
    _fmpz_vec_sub(v, v, res + 2 * m, f1len + g1len - 1);
    _fmpz_vec_add(res + m, res + m, v, vlen);

    _fmpz_vec_clear(t, tlen + (squaring ? 0 : ulen) + vlen);
}

/* Mulders' splitting point for short products of length n (measured to
   be flat between 0.6 and 0.8 here, with full products 10-25% slower) */
#define MULDERS_SPLIT(n) (((n) * 7) / 10)

/*
    Low product: coefficients [0, n) of f g, assuming flen, glen <= n.
    A full product of the low parts f0 g0 (split at k > n/2), plus
    recursive low products of the cross terms (which coincide when
    squaring).
*/
static void
_mullow(fmpz * res, const fmpz * f, slong flen, const fmpz * g, slong glen, slong n)
{
    slong L = flen + glen - 1;
    slong i, k, k1, k2, m, tlen;
    fmpz * t;
    int squaring = (f == g && flen == glen);

    FLINT_ASSERT(flen <= n && glen <= n && flen >= 1 && glen >= 1);

    if (L <= n)
    {
        _kar_mul(res, f, flen, g, glen);
        _fmpz_vec_zero(res + L, n - L);
        return;
    }

    if (_fmpz_poly_mulmid_toom_scalar_nmax(res, f, flen, g, glen, 0, n, _toom_scalar_nmax2(f, flen, g, glen)))
        return;

    /* nearly full: compute the full product */
    if (L - n <= n / 8)
    {
        t = _fmpz_vec_init(L);
        _kar_mul(t, f, flen, g, glen);
        for (i = 0; i < n; i++)
            fmpz_swap(res + i, t + i);
        _fmpz_vec_clear(t, L);
        return;
    }

    k = MULDERS_SPLIT(n);
    k = FLINT_MAX(k, (n + 1) / 2);
    k1 = FLINT_MIN(k, flen);
    k2 = FLINT_MIN(k, glen);
    m = n - k;

    /* f0 g0, keeping n coefficients */
    tlen = k1 + k2 - 1;
    if (tlen <= n)
    {
        _kar_mul(res, f, k1, g, k2);
        _fmpz_vec_zero(res + tlen, n - tlen);
    }
    else
    {
        t = _fmpz_vec_init(tlen);
        _kar_mul(t, f, k1, g, k2);
        for (i = 0; i < n; i++)
            fmpz_swap(res + i, t + i);
        _fmpz_vec_clear(t, tlen);
    }

    t = _fmpz_vec_init(m);

    /* x^k f1 g0, low part: g0 truncated to m coefficients */
    if (flen > k)
    {
        _mullow(t, f + k, flen - k, g, FLINT_MIN(glen, m), m);
        _fmpz_vec_add(res + k, res + k, t, m);
        /* when squaring, the other cross term is the same */
        if (squaring)
            _fmpz_vec_add(res + k, res + k, t, m);
    }

    /* x^k f0 g1, low part */
    if (glen > k && !squaring)
    {
        _mullow(t, f, FLINT_MIN(flen, m), g + k, glen - k, m);
        _fmpz_vec_add(res + k, res + k, t, m);
    }

    _fmpz_vec_clear(t, m);
}

/* High product: coefficients [nlo, flen + glen - 1) of f g, via reversal.
   The reversed operands are shallow (handle) copies, used read-only. */
static void
_mulhigh(fmpz * res, const fmpz * f, slong flen, const fmpz * g, slong glen, slong nlo)
{
    slong L = flen + glen - 1;
    slong n = L - nlo;
    slong i;
    fmpz * fr, * gr, * t;
    int squaring = (f == g && flen == glen);
    TMP_INIT;

    TMP_START;
    fr = TMP_ALLOC(sizeof(fmpz) * (flen + (squaring ? 0 : glen)));

    for (i = 0; i < flen; i++)
        fr[i] = f[flen - 1 - i];

    if (squaring)
    {
        gr = fr;
    }
    else
    {
        gr = fr + flen;
        for (i = 0; i < glen; i++)
            gr[i] = g[glen - 1 - i];
    }

    t = _fmpz_vec_init(n);
    _mullow(t, fr, FLINT_MIN(flen, n), gr, FLINT_MIN(glen, n), n);

    for (i = 0; i < n; i++)
        fmpz_swap(res + i, t + n - 1 - i);

    _fmpz_vec_clear(t, n);
    TMP_END;
}

/*
    Karatsuba middle product (Hanrot, Quercia, Zimmermann): sets res to
    the n coefficients [n - 1, 2n - 1) of a b where a has 2n - 1 and b
    has n coefficients.
*/
static void
_mulmid_kar(fmpz * res, const fmpz * a, const fmpz * b, slong n)
{
    slong i, k;
    fmpz * t1, * t2, * t3, * p1;

    if (_fmpz_poly_mulmid_toom_scalar_nmax(res, a, 2 * n - 1, b, n, n - 1, 2 * n - 1, _toom_scalar_nmax2(a, 2 * n - 1, b, n)))
        return;

    if (n % 2 == 1)
    {
        /* peel off the last coefficient of b and the last output */
        _mulmid_kar(res, a + 1, b, n - 1);
        for (i = 0; i < n - 1; i++)
            fmpz_addmul(res + i, a + i, b + n - 1);
        fmpz_mul(res + n - 1, a + 2 * n - 2, b + 0);
        for (i = 1; i < n; i++)
            fmpz_addmul(res + n - 1, a + 2 * n - 2 - i, b + i);
        return;
    }

    k = n / 2;

    t1 = _fmpz_vec_init(k + 2 * (2 * k - 1) + k);
    t2 = t1 + k;
    t3 = t2 + (2 * k - 1);
    p1 = t3 + (2 * k - 1);

    _fmpz_vec_add(t1, b, b + k, k);                     /* b0 + b1 */
    _fmpz_vec_sub(t2, a, a + k, 2 * k - 1);             /* a0 - a1 */
    _fmpz_vec_sub(t3, a + 2 * k, a + k, 2 * k - 1);     /* a2 - a1 */

    /* P1 = MP(a1, b0 + b1) */
    _mulmid_kar(p1, a + k, t1, k);
    /* res_lo = P1 + MP(a0 - a1, b1) */
    _mulmid_kar(res, t2, b + k, k);
    _fmpz_vec_add(res, res, p1, k);
    /* res_hi = P1 + MP(a2 - a1, b0) */
    _mulmid_kar(res + k, t3, b, k);
    _fmpz_vec_add(res + k, res + k, p1, k);

    _fmpz_vec_clear(t1, k + 2 * (2 * k - 1) + k);
}

void
_fmpz_poly_mulmid_toom_karatsuba(fmpz * res, const fmpz * poly1, slong len1,
                                 const fmpz * poly2, slong len2,
                                 slong nlo, slong nhi)
{
    const fmpz * f = poly1, * g = poly2;
    slong flen = len1, glen = len2;
    slong L, len, n, mlo, mhi, i;

    FLINT_ASSERT(flen >= 1 && glen >= 1);
    FLINT_ASSERT(0 <= nlo && nlo < nhi && nhi <= flen + glen - 1);

    /* trim the inputs to the coefficients that contribute to the window */
    flen = FLINT_MIN(flen, nhi);
    glen = FLINT_MIN(glen, nhi);
    {
        slong nlo2 = (flen + glen - 1) - nlo;

        if (flen > nlo2)
        {
            slong trunc = flen - nlo2;
            f += trunc;
            flen -= trunc;
            nlo -= trunc;
            nhi -= trunc;
        }

        if (glen > nlo2)
        {
            slong trunc = glen - nlo2;
            g += trunc;
            glen -= trunc;
            nlo -= trunc;
            nhi -= trunc;
        }
    }

    L = flen + glen - 1;
    len = nhi - nlo;

    if (flen < glen)
    {
        FLINT_SWAP(const fmpz *, f, g);
        FLINT_SWAP(slong, flen, glen);
    }

    if (nlo == 0 && nhi == L)
    {
        _kar_mul(res, f, flen, g, glen);
    }
    else if (_fmpz_poly_mulmid_toom_scalar_nmax(res, f, flen, g, glen, nlo, nhi, _toom_scalar_nmax2(f, flen, g, glen)))
    {
    }
    else if (nlo == 0)
    {
        _mullow(res, f, flen, g, glen, nhi);
    }
    else if (nhi == L)
    {
        _mulhigh(res, f, flen, g, glen, nlo);
    }
    else
    {
        /* Middle product region: coefficients [n - 1, flen) of f g
           (n = glen) involve all n coefficients of g and can be computed
           in blocks of n by the Karatsuba middle product. The edges are
           handled recursively (they trim to strictly smaller shapes). */
        n = glen;
        mlo = FLINT_MAX(nlo, n - 1);
        mhi = FLINT_MIN(nhi, flen);

        if (mhi - n < n - 1 || mhi - mlo < n / 2)
        {
            /* not a middle product shape: full product or classical */
            if (len >= L / 4)
            {
                fmpz * t = _fmpz_vec_init(L);
                _kar_mul(t, f, flen, g, glen);
                for (i = 0; i < len; i++)
                    fmpz_swap(res + i, t + nlo + i);
                _fmpz_vec_clear(t, L);
            }
            else
            {
                _fmpz_poly_mulmid_classical(res, f, flen, g, glen, nlo, nhi);
            }
        }
        else
        {
            if (nlo < mlo)
                _fmpz_poly_mulmid_toom_karatsuba(res, f, flen, g, glen, nlo, mlo);
            if (nhi > mhi)
                _fmpz_poly_mulmid_toom_karatsuba(res + (mhi - nlo), f, flen, g, glen, mhi, nhi);

            for (i = mlo; i < mhi; i += n)
            {
                if (mhi - i >= n)
                {
                    _mulmid_kar(res + (i - nlo), f + (i - (n - 1)), g, n);
                }
                else
                {
                    /* last partial block: compute a full block ending at mhi */
                    slong w = mhi - i;
                    slong start = mhi - n;
                    slong j;
                    fmpz * t = _fmpz_vec_init(n);
                    _mulmid_kar(t, f + (start - (n - 1)), g, n);
                    for (j = 0; j < w; j++)
                        fmpz_swap(res + (i - nlo) + j, t + (n - w) + j);
                    _fmpz_vec_clear(t, n);
                }
            }
        }
    }
}

void
_fmpz_poly_mul_toom_karatsuba(fmpz * res, const fmpz * poly1, slong len1,
                              const fmpz * poly2, slong len2)
{
    _kar_mul(res, poly1, len1, poly2, len2);
}

void
fmpz_poly_mul_toom_karatsuba(fmpz_poly_t res, const fmpz_poly_t poly1,
                             const fmpz_poly_t poly2)
{
    slong len1 = poly1->length, len2 = poly2->length;
    slong N = len1 + len2 - 1;

    if (len1 == 0 || len2 == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, N);
        _fmpz_poly_mul_toom_karatsuba(t->coeffs, poly1->coeffs, len1,
                                      poly2->coeffs, len2);
        _fmpz_poly_set_length(t, N);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
        return;
    }

    fmpz_poly_fit_length(res, N);
    _fmpz_poly_mul_toom_karatsuba(res->coeffs, poly1->coeffs, len1,
                                  poly2->coeffs, len2);
    _fmpz_poly_set_length(res, N);
}

void
fmpz_poly_mulmid_toom_karatsuba(fmpz_poly_t res, const fmpz_poly_t poly1,
                                const fmpz_poly_t poly2, slong nlo, slong nhi)
{
    slong len1 = poly1->length, len2 = poly2->length;
    slong w;

    FLINT_ASSERT(nlo >= 0);
    FLINT_ASSERT(nhi >= 0);

    if (len1 == 0 || len2 == 0 || nlo >= FLINT_MIN(nhi, len1 + len2 - 1))
    {
        fmpz_poly_zero(res);
        return;
    }

    nhi = FLINT_MIN(nhi, len1 + len2 - 1);
    w = nhi - nlo;

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t t;

        fmpz_poly_init2(t, w);
        _fmpz_poly_mulmid_toom_karatsuba(t->coeffs, poly1->coeffs, len1,
                                         poly2->coeffs, len2, nlo, nhi);
        _fmpz_poly_set_length(t, w);
        _fmpz_poly_normalise(t);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
        return;
    }

    fmpz_poly_fit_length(res, w);
    _fmpz_poly_mulmid_toom_karatsuba(res->coeffs, poly1->coeffs, len1,
                                     poly2->coeffs, len2, nlo, nhi);
    _fmpz_poly_set_length(res, w);
    _fmpz_poly_normalise(res);
}
