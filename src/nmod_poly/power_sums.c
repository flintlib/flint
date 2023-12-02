/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_poly.h"

void
_nmod_poly_power_sums(mp_ptr res, mp_srcptr poly, slong len, slong n,
                      nmod_t mod)
{
    if (10 * n >= len + 75)
        _nmod_poly_power_sums_schoenhage(res, poly, len, n, mod);
    else
        _nmod_poly_power_sums_naive(res, poly, len, n, mod);
}

void
nmod_poly_power_sums(nmod_poly_t res, const nmod_poly_t poly, slong n)
{
    slong len = poly->length;
    size_t i = 0;

    if (len == 0)
    {
        flint_throw(FLINT_ERROR, "(nmod_poly_power_sums_naive): Zero polynomial.\n");
    }
    while (poly->coeffs[i] == 0)
        i++;

    if (n <= 0 || len == 1)
        nmod_poly_zero(res);
    else if (len == i + 1)
    {
        nmod_poly_fit_length(res, 1);
        _nmod_poly_set_length(res, 1);
        NMOD_RED(res->coeffs[0], len - 1, poly->mod);
    }
    else
    {
        if (*nmod_poly_lead(poly) != 1)
        {
            nmod_poly_t t;
            nmod_poly_init_preinv(t, poly->mod.n, poly->mod.ninv);
            nmod_poly_make_monic(t, poly);
            nmod_poly_fit_length(res, n);
            _nmod_poly_power_sums(res->coeffs, t->coeffs + i,
                                  len - i, n, t->mod);
            nmod_poly_clear(t);
        }
        else if (poly == res)
        {
            nmod_poly_t t;
            nmod_poly_init_preinv(t, poly->mod.n, poly->mod.ninv);
            nmod_poly_fit_length(t, n);
            _nmod_poly_power_sums(t->coeffs, poly->coeffs + i,
                                  len - i, n, t->mod);
            nmod_poly_swap(t, res);
            nmod_poly_clear(t);
        }
        else
        {
            nmod_poly_fit_length(res, n);
            _nmod_poly_power_sums(res->coeffs, poly->coeffs + i,
                                  len - i, n, poly->mod);
        }
        if (i)
            NMOD_RED(res->coeffs[0], len - 1, poly->mod);
        _nmod_poly_set_length(res, n);
        _nmod_poly_normalise(res);
    }
}
