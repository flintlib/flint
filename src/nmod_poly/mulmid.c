/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2021, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "impl.h"

/*
    After trimming, the middle product computes nhi - nlo coefficients from
    operands of lengths len1 >= len2. It is the transpose of the product of
    lengths len2 and nhi - nlo, and costs about the same as that product
    with each algorithm except KS, which computes the full len1 x len2
    product. So fft_small is used where _nmod_poly_mul would use it for
    lengths len2 and nhi - nlo, and also from shorter lengths for some
    moduli where mul has KS variants in between; otherwise classical is
    used, or KS for small moduli and outputs much longer than len2.
*/
int
_nmod_poly_mulmid_want_fft_small(slong len2, slong m, nmod_t mod)
{
#if FLINT_HAVE_FFT_SMALL
    slong s = FLINT_MIN(len2, m);
    slong l = FLINT_MAX(len2, m);
    slong bits = NMOD_BITS(mod);

    /* fft_mul_tab is the crossover of mul with KS, KS2 and KS4; with
       classical as the alternative, fft_small is better from about 100
       for 9-15 bits and from 32 bits (22-31 bits: the table is kept) */
    return _nmod_poly_mullow_want_fft_small(l, s, l + s - 1, 0, mod)
            || (((bits >= 9 && bits <= 15) || bits >= 32)
                && FLINT_MIN(l, 2 * s) >= 100);
#else
    (void) len2;
    (void) m;
    (void) mod;
    return 0;
#endif
}

/* the dispatch after the IFMA kernel and the Karatsuba (see impl.h) */
FLINT_STATIC_NOINLINE void
_nmod_poly_mulmid_other(nn_ptr res, nn_srcptr poly1, slong len1,
                            nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
{
    slong m;
    int rows;

    /* trim, longer operand first (fft_small reads the second length) */
    _nmod_poly_mulmid_trim(&poly1, &len1, &poly2, &len2, &nlo, &nhi);

    m = nhi - nlo;

    /* rows (see impl.h), with twice the ratio of full and low products:
       here every dot product has len2 terms, while theirs are shorter near
       the ends; full products go to _nmod_poly_mul, which has its own rule */
    rows = m >= NMOD_VEC_SIMD_MIN_LEN && !(nlo == 0 && nhi == len1 + len2 - 1)
            && _nmod_poly_want_rows(_nmod_poly_mulmid_rows_min_len(len1, len2, nlo, nhi),
                                    len2, 2 * _NMOD_POLY_ROWS_RATIO(mod), mod)
            && !(poly1 == poly2 && len1 == len2);

    if (len2 <= 5 || m <= 5)
    {
        /* with NEON the rows in one word come first: on the M4 they beat
           the rows by 1.3-1.7x for len2 = 3 to 5 (Nov 2026) */
        if (rows && !(_NMOD_POLY_NEON && _nmod_poly_mulmid_want_one_word(len2, m, mod)))
            _nmod_poly_mulmid_rows(res, poly1, len1, poly2, len2, nlo, nhi, mod);
        else if (_nmod_poly_mulmid_want_one_word(len2, m, mod)
                 && !(nlo == 0 && nhi == len1 + len2 - 1)
                 && !(poly1 == poly2 && len1 == len2))
            _nmod_poly_mulmid_rows_one_word(res, poly1, len1, poly2, len2, nlo, nhi, mod);
        else
            _nmod_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi, mod);
        return;
    }

    if (nlo == 0 && nhi == len1 + len2 - 1)
    {
        _nmod_poly_mul(res, poly1, len1, poly2, len2, mod);
        return;
    }

    /* rows in one word (see impl.h), in place of the rows, the classical
       code, KS and fft_small */
    if (_nmod_poly_mulmid_want_one_word(len2, m, mod)
            && !(poly1 == poly2 && len1 == len2))
    {
        _nmod_poly_mulmid_rows_one_word(res, poly1, len1, poly2, len2, nlo, nhi, mod);
        return;
    }

#if FLINT_HAVE_FFT_SMALL
    if (_nmod_poly_mulmid_want_fft_small(len2, m, mod))
    {
        _nmod_poly_mulmid_fft_small(res, poly1, len1, poly2, len2, nlo, nhi, mod);
        return;
    }

    slong bits = NMOD_BITS(mod);

    /* KS computes the whole product: worth it for tiny moduli, for 16
       bits or less when the output is much longer than len2, and with NEON
       up to 12 bits for len2 > 10 or squares from 8 outputs (on the M4 KS
       beats the classical code there, which leads for other products with
       len2 <= 10 below 16 outputs, Nov 2026) */
    if (len2 >= 8 && ((bits <= 8 && m >= len2) || (bits <= 16 && m >= 8 * len2)
                      || (_NMOD_POLY_NEON && bits <= 12
                          && (len2 > 10 || (poly1 == poly2 && len1 == len2))
                          && m >= FLINT_MAX(8, (len2 + 3) / 4))))
        _nmod_poly_mulmid_KS(res, poly1, len1, poly2, len2, nlo, nhi, mod);
    else if (rows)
        _nmod_poly_mulmid_rows(res, poly1, len1, poly2, len2, nlo, nhi, mod);
    else
        _nmod_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi, mod);
#else
    {
        slong bits = NMOD_BITS(mod);
        slong n = FLINT_MIN(nhi, len1 + len2 - 1 - nlo);

        if (rows)
            _nmod_poly_mulmid_rows(res, poly1, len1, poly2, len2, nlo, nhi, mod);
        else if (n < 10 + bits * bits / 10)
            _nmod_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi, mod);
        else
            _nmod_poly_mulmid_KS(res, poly1, len1, poly2, len2, nlo, nhi, mod);
    }
#endif
}

/* as in _nmod_poly_mul, for the trimmed operands */
FLINT_STATIC_NOINLINE void
_nmod_poly_mulmid_kernels_n(nn_ptr res, nn_srcptr poly1, slong len1,
                            nn_srcptr poly2, slong len2, slong nlo, slong nhi,
                            nmod_t mod)
{
    int kind = _nmod_poly_product_kind(poly1, len1, poly2, len2, nlo, nhi);

    switch (_nmod_poly_dispatch_n(len1, len2, nhi - nlo, kind, mod))
    {
        case _NMOD_POLY_DISPATCH_KARATSUBA:
            _nmod_poly_mulmid_karatsuba(res, poly1, len1, poly2, len2, nlo,
                                        nhi, mod);
            break;
        case _NMOD_POLY_DISPATCH_IFMA:
            _nmod_poly_mulmid_ifma_ptr(res, poly1, len1, poly2, len2, nlo,
                                       nhi, &mod);
            break;
        default:
            _nmod_poly_mulmid_other(res, poly1, len1, poly2, len2, nlo, nhi,
                                    mod);
    }
}

FLINT_STATIC_NOINLINE void
_nmod_poly_mulmid_kernels(nn_ptr res, nn_srcptr poly1, slong len1,
                          nn_srcptr poly2, slong len2, slong nlo, slong nhi,
                          nmod_t mod)
{
    int kind = _nmod_poly_product_kind(poly1, len1, poly2, len2, nlo, nhi);

    switch (_nmod_poly_dispatch(len1, len2, nhi - nlo, kind, NMOD_BITS(mod)))
    {
        case _NMOD_POLY_DISPATCH_IFMA:
            _nmod_poly_mulmid_ifma(res, poly1, len1, poly2, len2, nlo, nhi,
                                   mod);
            break;
        case _NMOD_POLY_DISPATCH_KARATSUBA:
            _nmod_poly_mulmid_karatsuba(res, poly1, len1, poly2, len2, nlo,
                                        nhi, mod);
            break;
        case _NMOD_POLY_DISPATCH_BY_N:
            _nmod_poly_mulmid_kernels_n(res, poly1, len1, poly2, len2, nlo,
                                        nhi, mod);
            break;
        default:
            _nmod_poly_mulmid_other(res, poly1, len1, poly2, len2, nlo, nhi,
                                    mod);
    }
}

void
_nmod_poly_mulmid(nn_ptr res, nn_srcptr poly1, slong len1,
                            nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
{
    /* short operands or few outputs, or no kernel in this build (see
       impl.h): at once the other algorithms (the Karatsuba from 12
       coefficients) */
#if NMOD_POLY_HAVE_IFMA
    if (nhi - nlo <= 5
            || FLINT_MIN(len1, len2) < _NMOD_POLY_IFMA_MIN_LEN(NMOD_BITS(mod)))
#else
    if (!_NMOD_POLY_HAVE_KERNELS
            || FLINT_MIN(FLINT_MIN(len1, len2), nhi - nlo) < 12)
#endif
    {
        _nmod_poly_mulmid_other(res, poly1, len1, poly2, len2, nlo, nhi,
                                mod);
        return;
    }

    _nmod_poly_mulmid_trim(&poly1, &len1, &poly2, &len2, &nlo, &nhi);

    /* as in _nmod_poly_mul (_nmod_poly_mulmid_other trims the operands
       again, which changes nothing) */
    _nmod_poly_mulmid_kernels(res, poly1, len1, poly2, len2, nlo, nhi, mod);
}

void
nmod_poly_mulmid(nmod_poly_t res,
                           const nmod_poly_t poly1, const nmod_poly_t poly2,
                           slong nlo, slong nhi)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len;

    FLINT_ASSERT(nlo >= 0);
    FLINT_ASSERT(nhi >= 0);

    if (len1 == 0 || len2 == 0 || nlo >= FLINT_MIN(nhi, len1 + len2 - 1))
    {
        nmod_poly_zero(res);
        return;
    }

    nhi = FLINT_MIN(nhi, len1 + len2 - 1);
    len = nhi - nlo;

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;
        nmod_poly_init2_preinv(temp, poly1->mod.n, poly1->mod.ninv, len);
        _nmod_poly_mulmid(temp->coeffs, poly1->coeffs,
                                    len1, poly2->coeffs, len2, nlo, nhi, poly1->mod);
        nmod_poly_swap(res, temp);
        nmod_poly_clear(temp);
    }
    else
    {
        nmod_poly_fit_length(res, len);
        _nmod_poly_mulmid(res->coeffs, poly1->coeffs,
                                    len1, poly2->coeffs, len2, nlo, nhi, poly1->mod);
    }

    res->length = len;
    _nmod_poly_normalise(res);
}
