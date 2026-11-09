/*
    Copyright (C) 2011, 2021 Fredrik Johansson

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

void _nmod_poly_integral(nn_ptr res, nn_srcptr poly, slong len, nmod_t mod)
{
    /* n = 1: every coefficient is 0 (and k - 1 would not be reduced below) */
    if (mod.n == 1)
    {
        _nmod_vec_zero(res, len);
        return;
    }

    if (len > 2)
    {
        slong k;
        ulong t, u;

        res[len - 1] = poly[len - 2];
        t = len - 1;
        for (k = len - 2; k >= 2; k--)
        {
            res[k] = n_mulmod2_preinv(poly[k - 1], t, mod.n, mod.ninv);
            umul_ppmm(u, t, t, k);
            if (u != 0)
                t = n_ll_mod_preinv(u, t, mod.n, mod.ninv);
        }

        if (t >= mod.n)
            t = n_mod2_preinv(t, mod.n, mod.ninv);
        t = n_invmod(t, mod.n);

        /* t is reduced, and k - 1 < len - 1 < n: a prime factor of n >= 2 is
           larger than len - 1, or the inverse fails */
        res[2] = nmod_mul(res[2], t, mod);
        t = n_addmod(t, t, mod.n);

        if (len >= 4)
        {
            res[3] = nmod_mul(res[3], t, mod);

            for (k = 4; k < len; k++)
            {
                t = nmod_mul(t, k - 1, mod);
                res[k] = nmod_mul(res[k], t, mod);
            }
        }
    }

    if (len >= 2)
        res[1] = poly[0];

    res[0] = 0;
}

void nmod_poly_integral(nmod_poly_t x_int, const nmod_poly_t x)
{
    nmod_poly_fit_length(x_int, x->length + 1);
    _nmod_poly_integral(x_int->coeffs, x->coeffs, x->length + 1, x->mod);
    x_int->length = x->length + 1;
	_nmod_poly_normalise(x_int);
}
