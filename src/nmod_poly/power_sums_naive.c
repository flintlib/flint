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
_nmod_poly_power_sums_naive(mp_ptr res, mp_srcptr poly, slong len, slong n,
                            nmod_t mod)
{
    slong i, k;

    NMOD_RED(res[0], len - 1, mod);
    for (k = 1; k < FLINT_MIN(n, len); k++)
    {
        res[k] = nmod_mul(poly[len - 1 - k], k, mod);
        for (i = 1; i < k; i++)
            res[k] =
                nmod_add(res[k], nmod_mul(poly[len - 1 - k + i], res[i], mod),
                         mod);
        res[k] = nmod_neg(res[k], mod);
    }
    for (k = len; k < n; k++)
    {
        res[k] = 0;
        for (i = k - len + 1; i < k; i++)
            res[k] =
                nmod_add(res[k], nmod_mul(poly[len - 1 - k + i], res[i], mod),
                         mod);
        res[k] = nmod_neg(res[k], mod);
    }
}

void
nmod_poly_power_sums_naive(nmod_poly_t res, const nmod_poly_t poly, slong n)
{
    if (poly->length == 0)
    {
        flint_throw(FLINT_ERROR, "(nmod_poly_power_sums_naive): Zero polynomial.\n");
    }
    else if (n <= 0 || poly->length == 1)
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
            _nmod_poly_power_sums_naive(res->coeffs, t->coeffs,
                                        t->length, n, t->mod);
            nmod_poly_clear(t);
        }
        else if (poly == res)
        {
            nmod_poly_t t;
            nmod_poly_init_preinv(t, poly->mod.n, poly->mod.ninv);
            nmod_poly_fit_length(t, n);
            _nmod_poly_power_sums_naive(t->coeffs, poly->coeffs,
                                        poly->length, n, t->mod);
            nmod_poly_swap(t, res);
            nmod_poly_clear(t);
        }
        else
        {
            nmod_poly_fit_length(res, n);
            _nmod_poly_power_sums_naive(res->coeffs, poly->coeffs,
                                        poly->length, n, poly->mod);
        }
        _nmod_poly_set_length(res, n);
        _nmod_poly_normalise(res);
    }
}
