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

void
_nmod_poly_mulmid(nn_ptr res, nn_srcptr poly1, slong len1,
                            nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
{
    slong bits;
    slong n;

    len1 = FLINT_MIN(len1, nhi);
    len2 = FLINT_MIN(len2, nhi);

    if (len1 <= 5 || len2 <= 5 || nhi - nlo <= 5)
    {
        _nmod_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi, mod);
        return;
    }

    if (nlo == 0 && nhi == len1 + len2 - 1)
    {
        if (len1 >= len2)
            _nmod_poly_mul(res, poly1, len1, poly2, len2, mod);
        else
            _nmod_poly_mul(res, poly2, len2, poly1, len1, mod);
        return;
    }

    n = FLINT_MIN(nhi, len1 + len2 - 1 - nlo);

    if (_nmod_poly_mullow_want_fft_small(len1, len2, n, (poly1 == poly2 && len1 == len2), mod))
    {
        _nmod_poly_mulmid_fft_small(res, poly1, len1, poly2, len2, nlo, nhi, mod);
        return;
    }

    bits = NMOD_BITS(mod);

    if (n < 10 + bits * bits / 10)
        _nmod_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi, mod);
    else
        _nmod_poly_mulmid_KS(res, poly1, len1, poly2, len2, nlo, nhi, mod);
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
