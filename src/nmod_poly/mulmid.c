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
void
_nmod_poly_mulmid(nn_ptr res, nn_srcptr poly1, slong len1,
                            nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
{
    slong t, m;
#if FLINT_HAVE_FFT_SMALL
    slong s, l, bits;
#endif

    len1 = FLINT_MIN(len1, nhi);
    len2 = FLINT_MIN(len2, nhi);

    /* drop the low coefficients that no output coefficient uses */
    t = nlo - (len2 - 1);
    if (t > 0)
    {
        poly1 += t;
        len1 -= t;
        nlo -= t;
        nhi -= t;
    }

    t = nlo - (len1 - 1);
    if (t > 0)
    {
        poly2 += t;
        len2 -= t;
        nlo -= t;
        nhi -= t;
    }

    /* longer operand first: fft_small reads the second length */
    if (len1 < len2)
    {
        FLINT_SWAP(nn_srcptr, poly1, poly2);
        FLINT_SWAP(slong, len1, len2);
    }

    m = nhi - nlo;

    if (len2 <= 5 || m <= 5)
    {
        _nmod_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi, mod);
        return;
    }

    if (nlo == 0 && nhi == len1 + len2 - 1)
    {
        _nmod_poly_mul(res, poly1, len1, poly2, len2, mod);
        return;
    }

#if FLINT_HAVE_FFT_SMALL
    s = FLINT_MIN(len2, m);
    l = FLINT_MAX(len2, m);
    bits = NMOD_BITS(mod);

    /* fft_mul_tab is the crossover of mul with KS, KS2 and KS4; with
       classical as the alternative, fft_small is better from about 100
       for 9-15 bits and from 32 bits (22-31 bits: the table is kept) */
    if (_nmod_poly_mullow_want_fft_small(l, s, l + s - 1, 0, mod)
            || (((bits >= 9 && bits <= 15) || bits >= 32)
                && FLINT_MIN(l, 2 * s) >= 100))
    {
        _nmod_poly_mulmid_fft_small(res, poly1, len1, poly2, len2, nlo, nhi, mod);
        return;
    }

    /* KS computes the whole product: worth it for tiny moduli, and for
       16 bits or less when the output is much longer than len2 */
    if (len2 >= 8 && ((bits <= 8 && m >= len2) || (bits <= 16 && m >= 8 * len2)))
        _nmod_poly_mulmid_KS(res, poly1, len1, poly2, len2, nlo, nhi, mod);
    else
        _nmod_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi, mod);
#else
    {
        slong bits = NMOD_BITS(mod);
        slong n = FLINT_MIN(nhi, len1 + len2 - 1 - nlo);

        if (n < 10 + bits * bits / 10)
            _nmod_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi, mod);
        else
            _nmod_poly_mulmid_KS(res, poly1, len1, poly2, len2, nlo, nhi, mod);
    }
#endif
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
