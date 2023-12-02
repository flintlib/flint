/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_compose_mod(mp_ptr res,
    mp_srcptr f, slong lenf, mp_srcptr g, mp_srcptr h, slong lenh, nmod_t mod)
{
    if (lenh < 8 || lenf >= lenh)
        _nmod_poly_compose_mod_horner(res, f, lenf, g, h, lenh, mod);
    else
        _nmod_poly_compose_mod_brent_kung(res, f, lenf, g, h, lenh, mod);
}

void
nmod_poly_compose_mod(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_poly_t poly2,
                    const nmod_poly_t poly3)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;

    mp_ptr ptr2;

    if (len3 == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_compose_mod). Division by zero.\n");
    }

    if (len1 == 0 || len3 == 1)
    {
        nmod_poly_zero(res);
        return;
    }

    if (len1 == 1)
    {
        nmod_poly_set(res, poly1);
        return;
    }

    if (res == poly3 || res == poly1)
    {
        nmod_poly_t tmp;
        nmod_poly_init_preinv(tmp, res->mod.n, res->mod.ninv);
        nmod_poly_compose_mod(tmp, poly1, poly2, poly3);
        nmod_poly_swap(tmp, res);
        nmod_poly_clear(tmp);
        return;
    }

    ptr2 = _nmod_vec_init(len);

    if (len2 <= len)
    {
        flint_mpn_copyi(ptr2, poly2->coeffs, len2);
        flint_mpn_zero(ptr2 + len2, len - len2);
    }
    else
    {
        _nmod_poly_rem(ptr2, poly2->coeffs, len2,
                             poly3->coeffs, len3, res->mod);
    }

    nmod_poly_fit_length(res, len);
    _nmod_poly_compose_mod(res->coeffs,
        poly1->coeffs, len1, ptr2, poly3->coeffs, len3, res->mod);
    res->length = len;
    _nmod_poly_normalise(res);

    _nmod_vec_clear(ptr2);
}
