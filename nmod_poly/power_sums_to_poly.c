/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_power_sums_to_poly(mp_ptr res, mp_srcptr poly, slong len,
                              nmod_t mod)
{
    if (mod.n <= 12 || poly[0] <= 10)
        _nmod_poly_power_sums_to_poly_naive(res, poly, len, mod);
    else
        _nmod_poly_power_sums_to_poly_schoenhage(res, poly, len, mod);
}

void
nmod_poly_power_sums_to_poly(nmod_poly_t res, const nmod_poly_t Q)
{
    if (Q->length == 0)
    {
        nmod_poly_fit_length(res, 1);
        res->coeffs[0] = 1;
        _nmod_poly_set_length(res, 1);
    }
    else
    {
        slong d = Q->coeffs[0];
        if (Q == res)
        {
            nmod_poly_t t;
            nmod_poly_init_preinv(t, Q->mod.n, Q->mod.ninv);
            nmod_poly_fit_length(t, d + 1);
            _nmod_poly_power_sums_to_poly(t->coeffs, Q->coeffs, Q->length,
                                          Q->mod);
            nmod_poly_swap(res, t);
            nmod_poly_clear(t);
        }
        else
        {
            nmod_poly_fit_length(res, d + 1);
            _nmod_poly_power_sums_to_poly(res->coeffs, Q->coeffs, Q->length,
                                          Q->mod);
        }
        _nmod_poly_set_length(res, d + 1);
        _nmod_poly_normalise(res);
    }
}
