/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "impl.h"

/* the dispatch after the IFMA kernel and the Karatsuba (see impl.h) */
FLINT_STATIC_NOINLINE void
_nmod_poly_mullow_other(nn_ptr res, nn_srcptr poly1, slong len1,
                             nn_srcptr poly2, slong len2, slong n, nmod_t mod)
{
    slong bits;
    int squaring, rows, one_word;

    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    squaring = (poly1 == poly2 && len1 == len2);

    /* rows of length len1 down to n - len2 + 1 (see impl.h); when the
       classical code accumulates in one word, only for len2 <= 4, as in
       _nmod_poly_mul */
    one_word = 2 * NMOD_BITS(mod) + FLINT_BIT_COUNT(len2) <= FLINT_BITS;
    rows = _nmod_poly_want_rows(n - len2 + 1, len2, _NMOD_POLY_ROWS_RATIO(mod), mod)
                && !squaring && (!one_word || len2 <= 4);

    if (len2 <= 5)
    {
        if (rows)
            _nmod_poly_mulmid_rows(res, poly1, len1, poly2, len2, 0, n, mod);
        else
            _nmod_poly_mullow_classical(res, poly1, len1, poly2, len2, n, mod);
        return;
    }

    if (n == len1 + len2 - 1)
    {
        _nmod_poly_mul(res, poly1, len1, poly2, len2, mod);
        return;
    }

    if (_nmod_poly_mullow_want_fft_small(len1, len2, n, squaring, mod))
    {
        _nmod_poly_mullow_fft_small(res, poly1, len1, poly2, len2, n, mod);
        return;
    }

    bits = NMOD_BITS(mod);

    /* the classical code below n = 10 + b^2/10, and when it accumulates
       in one word, from 16 bits at any length below fft_small (KS computes
       the whole product; 1.3-2.5x from 20 bits) and below for len2 < c(b)
       (see _nmod_poly_one_word_cutoff). With NEON, on the M4 KS wins from
       n = 8-14 below 24 bits, and for len2 from 6 below 16 bits, from 7
       above (Nov 2026): there the classical code only for n < 10, or
       len2 < 6 below 16 bits and len2 < 7 from 16 bits */
    if (rows)
        _nmod_poly_mulmid_rows(res, poly1, len1, poly2, len2, 0, n, mod);
    else if (_NMOD_POLY_NEON
             ? (n < (bits < 24 ? 10 : 10 + bits * bits / 10)
                || (bits < 24 && len2 < (bits < 16 ? 6 : 7)))
             : (n < 10 + bits * bits / 10
                || (one_word && (bits >= 16
                                 || len2 < _nmod_poly_one_word_cutoff(bits)))))
        _nmod_poly_mullow_classical(res, poly1, len1, poly2, len2, n, mod);
    else
        _nmod_poly_mullow_KS(res, poly1, len1, poly2, len2, n, mod);
}

/* as in _nmod_poly_mul, for len1, len2 <= n */
FLINT_STATIC_NOINLINE void
_nmod_poly_mullow_kernels_n(nn_ptr res, nn_srcptr poly1, slong len1,
                            nn_srcptr poly2, slong len2, slong n, nmod_t mod)
{
    int kind = _nmod_poly_product_kind(poly1, len1, poly2, len2, 0, n);

    switch (_nmod_poly_dispatch_n(len1, len2, n, kind, mod))
    {
        case _NMOD_POLY_DISPATCH_KARATSUBA:
            _nmod_poly_mulmid_karatsuba(res, poly1, len1, poly2, len2, 0, n, mod);
            break;
        case _NMOD_POLY_DISPATCH_IFMA:
            _nmod_poly_mulmid_ifma_ptr(res, poly1, len1, poly2, len2, 0, n, &mod);
            break;
        default:
            _nmod_poly_mullow_other(res, poly1, len1, poly2, len2, n, mod);
    }
}

FLINT_STATIC_NOINLINE void
_nmod_poly_mullow_kernels(nn_ptr res, nn_srcptr poly1, slong len1,
                          nn_srcptr poly2, slong len2, slong n, nmod_t mod)
{
    int kind = _nmod_poly_product_kind(poly1, len1, poly2, len2, 0, n);

    switch (_nmod_poly_dispatch(len1, len2, n, kind, NMOD_BITS(mod)))
    {
        case _NMOD_POLY_DISPATCH_IFMA:
            _nmod_poly_mulmid_ifma(res, poly1, len1, poly2, len2, 0, n, mod);
            break;
        case _NMOD_POLY_DISPATCH_KARATSUBA:
            _nmod_poly_mulmid_karatsuba(res, poly1, len1, poly2, len2, 0, n, mod);
            break;
        case _NMOD_POLY_DISPATCH_BY_N:
            _nmod_poly_mullow_kernels_n(res, poly1, len1, poly2, len2, n, mod);
            break;
        default:
            _nmod_poly_mullow_other(res, poly1, len1, poly2, len2, n, mod);
    }
}

void _nmod_poly_mullow(nn_ptr res, nn_srcptr poly1, slong len1,
                             nn_srcptr poly2, slong len2, slong n, nmod_t mod)
{
    slong s = FLINT_MIN(len2, n);

    /* as in _nmod_poly_mul */
    if (!_NMOD_POLY_HAVE_KERNELS
            || (s < _NMOD_POLY_DISPATCH_MIN_LEN(NMOD_BITS(mod))
                && !_nmod_poly_karatsuba_from_8(s, n,
                    _nmod_poly_product_kind(poly1, len1, poly2, len2, 0, n),
                    NMOD_BITS(mod))))
        _nmod_poly_mullow_other(res, poly1, len1, poly2, len2, n, mod);
    else
        _nmod_poly_mullow_kernels(res, poly1, FLINT_MIN(len1, n), poly2,
                                  FLINT_MIN(len2, n), n, mod);
}

void nmod_poly_mullow(nmod_poly_t res,
              const nmod_poly_t poly1, const nmod_poly_t poly2, slong trunc)
{
    slong len1, len2, len_out;

    len1 = poly1->length;
    len2 = poly2->length;

    len_out = poly1->length + poly2->length - 1;
    if (trunc > len_out)
        trunc = len_out;

    if (len1 == 0 || len2 == 0 || trunc == 0)
    {
        nmod_poly_zero(res);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;

        nmod_poly_init2(temp, poly1->mod.n, trunc);

        if (len1 >= len2)
            _nmod_poly_mullow(temp->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, trunc, poly1->mod);
        else
            _nmod_poly_mullow(temp->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, trunc, poly1->mod);

        nmod_poly_swap(temp, res);
        nmod_poly_clear(temp);
    } else
    {
        nmod_poly_fit_length(res, trunc);

        if (len1 >= len2)
            _nmod_poly_mullow(res->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, trunc, poly1->mod);
        else
            _nmod_poly_mullow(res->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, trunc, poly1->mod);
    }

    res->length = trunc;
    _nmod_poly_normalise(res);
}
