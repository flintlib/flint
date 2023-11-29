/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"

void
_nmod_poly_power_sums_schoenhage(mp_ptr res, mp_srcptr poly, slong len,
                                 slong n, nmod_t mod)
{
    mp_ptr a, b;

    a = (mp_ptr) flint_malloc((2 * len - 1) * sizeof(mp_limb_t));
    b = a + len;

    _nmod_poly_reverse(a, poly, len, len);
    _nmod_poly_derivative(b, poly, len, mod);
    _nmod_poly_reverse(b, b, len - 1, len - 1);

    _nmod_poly_div_series(res, b, len - 1, a, len, n, mod);

    flint_free(a);
}

void
nmod_poly_power_sums_schoenhage(nmod_poly_t res, const nmod_poly_t poly,
                                slong n)
{
    if (poly->length == 0)
    {
        flint_throw(FLINT_ERROR, "(nmod_poly_power_sums_schoenhage): Zero polynomial.\n");
    }
    else if ((n <= 0) || (poly->length == 1))
    {
        nmod_poly_zero(res);
    }
    else
    {
        if (*nmod_poly_lead(poly) != 1)
        {
            nmod_poly_t t;
            nmod_poly_init_preinv(t, poly->mod.n, poly->mod.ninv);
            nmod_poly_make_monic(t, poly);
            nmod_poly_fit_length(res, n);
            _nmod_poly_power_sums_schoenhage(res->coeffs, t->coeffs,
                                             t->length, n, t->mod);
            nmod_poly_clear(t);
        }
        else if (poly == res)
        {
            nmod_poly_t t;
            nmod_poly_init_preinv(t, poly->mod.n, poly->mod.ninv);
            nmod_poly_fit_length(t, n);
            _nmod_poly_power_sums_schoenhage(t->coeffs, poly->coeffs,
                                             poly->length, n, t->mod);
            nmod_poly_swap(t, res);
            nmod_poly_clear(t);
        }
        else
        {
            nmod_poly_fit_length(res, n);
            _nmod_poly_power_sums_schoenhage(res->coeffs, poly->coeffs,
                                             poly->length, n, poly->mod);
        }
        _nmod_poly_set_length(res, n);
        _nmod_poly_normalise(res);
    }
}
