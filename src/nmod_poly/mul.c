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
#include "nmod_poly.h"

#if FLINT_HAVE_FFT_SMALL

#include "fft_small.h"

/* todo: check unbalanced cutoffs */

static const short fft_mul_tab[] = {1326, 1326, 1095, 802, 674, 537, 330, 306, 290,
274, 200, 192, 182, 173, 163, 99, 97, 93, 90, 82, 80, 438, 414, 324, 393,
298, 298, 268, 187, 185, 176, 176, 168, 167, 158, 158, 97, 96, 93, 92, 89,
89, 85, 85, 80, 81, 177, 172, 163, 162, 164, 176, 171, 167, 167, 164, 163,
163, 160, 165, 95, 96, 90, 94, };

static const short fft_sqr_tab[] = {1420, 1420, 1353, 964, 689, 569, 407, 353, 321,
321, 292, 279, 200, 182, 182, 159, 159, 152, 145, 139, 723, 626, 626, 569,
597, 448, 542, 292, 292, 200, 191, 191, 182, 182, 166, 166, 166, 159, 159,
159, 152, 152, 145, 145, 93, 200, 191, 182, 182, 182, 182, 191, 191, 191,
182, 182, 174, 182, 182, 182, 152, 152, 152, 145, };

#endif

void _nmod_poly_mul(nn_ptr res, nn_srcptr poly1, slong len1,
                             nn_srcptr poly2, slong len2, nmod_t mod)
{
    slong bits, cutoff_len;

    if (len2 <= 5)
    {
        _nmod_poly_mul_classical(res, poly1, len1, poly2, len2, mod);
        return;
    }

    bits = NMOD_BITS(mod);
    cutoff_len = FLINT_MIN(len1, 2 * len2);

#if FLINT_HAVE_FFT_SMALL

    if (poly1 == poly2 && len1 == len2)
    {
        if (cutoff_len >= fft_sqr_tab[bits - 1])
        {
            _nmod_poly_mul_mid_default_mpn_ctx(res, 0, len1 + len2 - 1, poly1, len1, poly2, len2, mod);
            return;
        }
    }
    else
    {
        if (cutoff_len >= fft_mul_tab[bits - 1])
        {
            _nmod_poly_mul_mid_default_mpn_ctx(res, 0, len1 + len2 - 1, poly1, len1, poly2, len2, mod);
            return;
        }
    }

#endif

    if (3 * cutoff_len < 2 * FLINT_MAX(bits, 10))
        _nmod_poly_mul_classical(res, poly1, len1, poly2, len2, mod);
    else if (cutoff_len * bits < 800)
        _nmod_poly_mul_KS(res, poly1, len1, poly2, len2, 0, mod);
    else if (cutoff_len * (bits + 1) * (bits + 1) < 100000)
        _nmod_poly_mul_KS2(res, poly1, len1, poly2, len2, mod);
    else
        _nmod_poly_mul_KS4(res, poly1, len1, poly2, len2, mod);
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
