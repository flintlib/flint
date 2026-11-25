/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "impl.h"

/* whether _nmod_poly_mul_classical accumulates in one word (dot method
   _DOT1, not _DOT_POW2), up to the rounding of the bounds of
   _nmod_vec_dot_params for len2 <= 11 */
static int
_mul_classical_one_word(slong len2, nmod_t mod)
{
    return _nmod_poly_one_word(len2, mod) && (mod.n & (mod.n - 1)) != 0;
}

/*
    Whether to take the classical code, accumulating in one word, rather
    than KS (or KS2, KS4): for len1 >= 2 len2, len2 < c(b) (see
    _nmod_poly_one_word_cutoff); from 20 bits, also for balanced products
    with len2 < c(b) and squares with len2 < 8b - 96 (Cascade Lake, Nov
    2026: balanced up to 36-48 at 20 bits, 64-96 at 22-28 bits; squares
    64-80 at 20 bits, 112-160 at 24-28 bits); below 20 bits, KS beats the
    balanced classical code from length 6-12. One word allows at most 32
    bits. With NEON, only squares from 22 bits for len2 < 2b - 8: on the
    M4, KS is faster at the other lengths, and squares up to 32-56 from 22
    bits (at 20-21 bits the classical squares lead by at most 13% and lose
    in fq_nmod_pow, 0.83x in degree 16, Nov 2026).
*/
static int
_mul_want_classical_one_word(slong len1, slong len2, int squaring, nmod_t mod)
{
    slong bits = NMOD_BITS(mod), c;

    if (_NMOD_POLY_NEON)
        return squaring && bits >= 22 && len2 < 2 * bits - 8
            && _mul_classical_one_word(len2, mod);

    c = _nmod_poly_one_word_cutoff(bits);
    if (len1 < 2 * len2)
        c = (bits < 20) ? 0 : squaring ? 8 * bits - 96 : c;

    return len2 < c && _mul_classical_one_word(len2, mod);
}

/* the dispatch after the IFMA kernel and the Karatsuba (see impl.h) */
FLINT_STATIC_NOINLINE void
_nmod_poly_mul_other(nn_ptr res, nn_srcptr poly1, slong len1,
                             nn_srcptr poly2, slong len2, nmod_t mod)
{
    slong bits, cutoff_len;
    int squaring, rows;

    squaring = (poly1 == poly2 && len1 == len2);

    /* rows of length len1 (see impl.h); when the classical code
       accumulates in one word, only for len2 <= 4 (there 1.2-1.5x faster
       on Ice Lake, Zen 4 and the M4, 0.85-1.0x on a Cascade Lake, Nov
       2026); for powers of 2 below 2^16, KS beats the rows from len2 = 12
       (0.6-0.9x, Nov 2026), and with NEON KS or the classical code beat
       them from len2 = 7 up to 2^31 (on the M4, 0.6-0.9x, Nov 2026) */
    rows = _nmod_poly_want_rows(len1, len2, _NMOD_POLY_ROWS_RATIO(mod), mod)
                && !squaring && (len2 <= 4 || !_mul_classical_one_word(len2, mod))
                && (len2 <= 8 || NMOD_BITS(mod) > 16)
                && !(_NMOD_POLY_NEON && len2 > 6 && (mod.n & (mod.n - 1)) == 0);

    if (len2 <= 5)
    {
        if (rows)
            _nmod_poly_mulmid_rows(res, poly1, len1, poly2, len2, 0, len1 + len2 - 1, mod);
        else
            _nmod_poly_mul_classical(res, poly1, len1, poly2, len2, mod);
        return;
    }

    if (_nmod_poly_mullow_want_fft_small(len1, len2, len1 + len2 - 1, squaring, mod))
    {
        _nmod_poly_mullow_fft_small(res, poly1, len1, poly2, len2, len1 + len2 - 1, mod);
        return;
    }

    if (rows)
    {
        _nmod_poly_mulmid_rows(res, poly1, len1, poly2, len2, 0, len1 + len2 - 1, mod);
        return;
    }

    bits = NMOD_BITS(mod);
    cutoff_len = FLINT_MIN(len1, 2 * len2);

    if (3 * cutoff_len < 2 * FLINT_MAX(bits, 10)
            || _mul_want_classical_one_word(len1, len2, squaring, mod))
        _nmod_poly_mul_classical(res, poly1, len1, poly2, len2, mod);
    else if (cutoff_len * bits < 800)
        _nmod_poly_mul_KS(res, poly1, len1, poly2, len2, mod);
    else if (cutoff_len * (bits + 1) * (bits + 1) < 100000)
        _nmod_poly_mul_KS2(res, poly1, len1, poly2, len2, mod);
    else
        _nmod_poly_mul_KS4(res, poly1, len1, poly2, len2, mod);
}

/* where n decides (see impl.h): the IFMA kernel, with the one-word
   Karatsuba before it for longer operands, else the Karatsuba if its rule
   takes the product, else the other algorithms */
FLINT_STATIC_NOINLINE void
_nmod_poly_mul_kernels_n(nn_ptr res, nn_srcptr poly1, slong len1,
                         nn_srcptr poly2, slong len2, nmod_t mod)
{
    slong m = len1 + len2 - 1;
    int kind = _nmod_poly_product_kind(poly1, len1, poly2, len2, 0, m);

    switch (_nmod_poly_dispatch_n(len1, len2, m, kind, mod))
    {
        case _NMOD_POLY_DISPATCH_KARATSUBA:
            _nmod_poly_mulmid_karatsuba(res, poly1, len1, poly2, len2, 0, m, mod);
            break;
        case _NMOD_POLY_DISPATCH_IFMA:
            _nmod_poly_mulmid_ifma_ptr(res, poly1, len1, poly2, len2, 0, m, &mod);
            break;
        default:
            _nmod_poly_mul_other(res, poly1, len1, poly2, len2, mod);
    }
}

/* from the lengths and the bit size of n (see impl.h) */
FLINT_STATIC_NOINLINE void
_nmod_poly_mul_kernels(nn_ptr res, nn_srcptr poly1, slong len1,
                       nn_srcptr poly2, slong len2, nmod_t mod)
{
    slong m = len1 + len2 - 1;
    int kind = _nmod_poly_product_kind(poly1, len1, poly2, len2, 0, m);

    switch (_nmod_poly_dispatch(len1, len2, m, kind, NMOD_BITS(mod)))
    {
        case _NMOD_POLY_DISPATCH_IFMA:
            _nmod_poly_mulmid_ifma(res, poly1, len1, poly2, len2, 0, m, mod);
            break;
        case _NMOD_POLY_DISPATCH_KARATSUBA:
            _nmod_poly_mulmid_karatsuba(res, poly1, len1, poly2, len2, 0, m, mod);
            break;
        case _NMOD_POLY_DISPATCH_BY_N:
            _nmod_poly_mul_kernels_n(res, poly1, len1, poly2, len2, mod);
            break;
        default:
            _nmod_poly_mul_other(res, poly1, len1, poly2, len2, mod);
    }
}

void _nmod_poly_mul(nn_ptr res, nn_srcptr poly1, slong len1,
                             nn_srcptr poly2, slong len2, nmod_t mod)
{
    /* short operands (all operands in a build without kernels) at once to
       the other algorithms, longer ones to _nmod_poly_mul_kernels, which
       decides from the bit size of n; both calls come last, so that the
       compiler can make them jumps (tail calls, see impl.h) */
    if (!_NMOD_POLY_HAVE_KERNELS
            || (len2 < _NMOD_POLY_DISPATCH_MIN_LEN(NMOD_BITS(mod))
                && !_nmod_poly_karatsuba_from_8(len2, len1 + len2 - 1,
                    _nmod_poly_product_kind(poly1, len1, poly2, len2, 0,
                                            len1 + len2 - 1),
                    NMOD_BITS(mod))))
        _nmod_poly_mul_other(res, poly1, len1, poly2, len2, mod);
    else
        _nmod_poly_mul_kernels(res, poly1, len1, poly2, len2, mod);
}

void nmod_poly_mul(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2)
{
    slong len1, len2, len_out;

    len1 = poly1->length;
    len2 = poly2->length;

    if (len1 == 0 || len2 == 0)
    {
        nmod_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;

        nmod_poly_init2(temp, poly1->mod.n, len_out);

        if (len1 >= len2)
            _nmod_poly_mul(temp->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, poly1->mod);
        else
            _nmod_poly_mul(temp->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, poly1->mod);

        nmod_poly_swap(temp, res);
        nmod_poly_clear(temp);
    } else
    {
        nmod_poly_fit_length(res, len_out);

        if (len1 >= len2)
            _nmod_poly_mul(res->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, poly1->mod);
        else
            _nmod_poly_mul(res->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, poly1->mod);
    }

    res->length = len_out;
    _nmod_poly_normalise(res);
}
