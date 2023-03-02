/*
    Copyright (C) 2011, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void _nmod_poly_integral(mp_ptr res, mp_srcptr poly, slong len, nmod_t mod)
{
    if (len > 2)
    {
        slong k;
        mp_limb_t t, u;

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

        res[2] = n_mulmod2_preinv(res[2], t, mod.n, mod.ninv);
        t = n_addmod(t, t, mod.n);

        if (len >= 4)
        {
            res[3] = n_mulmod2_preinv(res[3], t, mod.n, mod.ninv);

            for (k = 4; k < len; k++)
            {
                t = n_mulmod2_preinv(t, k - 1, mod.n, mod.ninv);
                res[k] = n_mulmod2_preinv(res[k], t, mod.n, mod.ninv);
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
