/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/******************************************************************************

    Authored 2016 by Daniel S. Roche; US Government work in the public domain.

******************************************************************************/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_poly.h"

slong
_fmpz_mod_poly_minpoly_bm(fmpz * poly, const fmpz * seq, slong len, const fmpz_mod_ctx_t ctx)
{
    fmpz *buf, *curpoly, *prevpoly;
    slong curlen, prevlen;
    fmpz_t disc;
    slong i, m;

    buf = _fmpz_vec_init(len + 1);
    curpoly = poly;
    _fmpz_vec_zero(curpoly, len + 1);
    prevpoly = buf;
    fmpz_init(disc);

    fmpz_one(curpoly + 0);
    curlen = 1;
    fmpz_one(prevpoly + 0);
    prevlen = 1;
    m = -1; /* m = last switching point */

    for (i = 0; i < len; ++i)
    {
        /* compute next discrepancy */
        _fmpz_vec_dot(disc, curpoly, seq + (i - curlen + 1), curlen);
        fmpz_mod_set_fmpz(disc, disc, ctx);

        if (fmpz_is_zero(disc))
        {
        }
        else if (i - m <= curlen - prevlen)
        {
            /* quick update; no switch, curlen doesn't change. */
            slong pos = (curlen - prevlen) - (i - m);
            _fmpz_vec_scalar_addmul_fmpz(curpoly + pos, prevpoly, prevlen, disc);
            _fmpz_mod_vec_set_fmpz_vec(curpoly + pos, curpoly + pos, prevlen, ctx);
        }
        else
        {
            /* switching update */
            slong pos = (i - m) - (curlen - prevlen);

            _fmpz_mod_vec_scalar_mul_fmpz_mod(prevpoly, prevpoly, prevlen, disc, ctx);
            _fmpz_mod_poly_add(prevpoly + pos, prevpoly + pos, FLINT_MAX(0, prevlen - pos), curpoly, curlen, ctx);
            prevlen = curlen + pos;

            fmpz_mod_neg(disc, disc, ctx);
            fmpz_mod_inv(disc, disc, ctx);
            _fmpz_mod_vec_scalar_mul_fmpz_mod(curpoly, curpoly, curlen, disc, ctx);

            FMPZ_VEC_SWAP(curpoly, curlen, prevpoly, prevlen);
            m = i;
        }
    }

    /* make curpoly monic and copy to poly if necessary */
    fmpz_mod_inv(disc, curpoly + (curlen - 1), ctx);
    _fmpz_mod_poly_scalar_mul_fmpz(poly, curpoly, curlen, disc, ctx);

    _fmpz_vec_clear(buf, len + 1);
    fmpz_clear(disc);

    return curlen;
}

void
fmpz_mod_poly_minpoly_bm(fmpz_mod_poly_t poly, const fmpz * seq, slong len,
                                                      const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_fit_length(poly, len + 1, ctx);
    poly->length = _fmpz_mod_poly_minpoly_bm(poly->coeffs, seq, len, ctx);
}

slong
_fmpz_mod_poly_minpoly_hgcd(fmpz * poly, const fmpz * seq, slong len, const fmpz_mod_ctx_t ctx)
{
    fmpz *buf, *f, *g, *A, *B;
    fmpz* M[4];
    slong lenM[4];
    slong buflen, lenA, lenB, len_poly, leng;
    int i;

    M[0] = poly;
    buflen = 7 * len + 5; /* for f, g, A, B, M[1], M[2], M[3] */
    buf = _fmpz_vec_init(buflen);
    f = buf;
    g = f + (len + 1);
    A = g + len;
    B = A + (len + 1);
    M[1] = B + len;
    M[2] = M[1] + (len + 1);
    M[3] = M[2] + (len + 1);

    /* f = x^len */
    fmpz_one(f + len);
    /* g = reversal of seq */
    for (i = 0; i < len; ++i) fmpz_set(g + i, seq + (len - i - 1));
    leng = len;
    FMPZ_VEC_NORM(g, leng);

    /* leng is invalid input for hgcd. todo: change hgcd to allow this? */
    if (leng == 0)
    {
        fmpz_one(M[0]);
        fmpz_one(M[3]);
        lenM[0] = lenM[3] = 1;
        lenM[1] = lenM[2] = 0;
        lenA = len + 1;
        _fmpz_vec_set(A, f, lenA);
        lenB = leng;
        _fmpz_vec_set(B, g, leng);
    }
    else
    {
        _fmpz_mod_poly_hgcd(M, lenM, A, &lenA, B, &lenB, f, len + 1, g, leng, ctx);
    }

    len_poly = lenM[0];

    /* one more step may be necessary */
    if (len_poly <= lenB)
    {
        slong quo_len = lenA - lenB + 1;
        fmpz_mod_inv(buf, B + (lenB - 1), ctx);
        _fmpz_mod_poly_divrem(M[2], M[3], A, lenA, B, lenB, buf, ctx);

        if (len_poly >= quo_len)
        {
            _fmpz_mod_poly_mul(M[3], poly, len_poly, M[2], quo_len, ctx);
        }
        else
        {
            _fmpz_mod_poly_mul(M[3], M[2], quo_len, poly, len_poly, ctx);
        }

        len_poly += quo_len - 1;
        _fmpz_mod_poly_add(poly, M[3], len_poly, M[1], lenM[1], ctx);
    }

    /* make poly monic */
    fmpz_mod_inv(buf, poly + (len_poly - 1), ctx);
    _fmpz_mod_poly_scalar_mul_fmpz(poly, poly, len_poly, buf, ctx);

    _fmpz_vec_clear(buf, buflen);

    return len_poly;
}

void
fmpz_mod_poly_minpoly_hgcd(fmpz_mod_poly_t poly, const fmpz * seq, slong len,
                                                      const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_fit_length(poly, len + 1, ctx);
    poly->length = _fmpz_mod_poly_minpoly_hgcd(poly->coeffs, seq, len, ctx);
}

slong
_fmpz_mod_poly_minpoly(fmpz * poly, const fmpz * seq, slong len, const fmpz_mod_ctx_t ctx)
{
    if (len < FLINT_MAX(200, 530-22*fmpz_size(fmpz_mod_ctx_modulus(ctx))))
        return _fmpz_mod_poly_minpoly_bm(poly, seq, len, ctx);
    else
        return _fmpz_mod_poly_minpoly_hgcd(poly, seq, len, ctx);
}

void
fmpz_mod_poly_minpoly(fmpz_mod_poly_t poly, const fmpz * seq, slong len,
                                                      const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_fit_length(poly, len+1, ctx);
    poly->length = _fmpz_mod_poly_minpoly(poly->coeffs, seq, len, ctx);
}
