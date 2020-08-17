/*
    Copyright (C) 2011 Fredrik Johansson

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


void
_nmod_poly_interpolate_nmod_vec_barycentric(mp_ptr poly,
                            mp_srcptr xs, mp_srcptr ys, slong n, nmod_t mod)
{
    mp_ptr P, Q, w;
    slong i, j;

    if (n == 1)
    {
        poly[0] = ys[0];
        return;
    }

    P = _nmod_vec_init(n + 1);
    Q = _nmod_vec_init(n);
    w = _nmod_vec_init(n);

    _nmod_poly_product_roots_nmod_vec(P, xs, n, mod);

    for (i = 0; i < n; i++)
    {
        w[i] = UWORD(1);
        for (j = 0; j < n; j++)
        {
            if (i != j)
                w[i] = nmod_mul(w[i], nmod_sub(xs[i], xs[j], mod), mod);
        }
        w[i] = n_invmod(w[i], mod.n);
    }

    _nmod_vec_zero(poly, n);

    for (i = 0; i < n; i++)
    {
        _nmod_poly_div_root(Q, P, n + 1, xs[i], mod);
        _nmod_vec_scalar_addmul_nmod(poly, Q, n,
            nmod_mul(w[i], ys[i], mod), mod);
    }

    _nmod_vec_clear(P);
    _nmod_vec_clear(Q);
    _nmod_vec_clear(w);
}

void
nmod_poly_interpolate_nmod_vec_barycentric(nmod_poly_t poly,
                                    mp_srcptr xs, mp_srcptr ys, slong n)
{
    if (n == 0)
    {
        nmod_poly_zero(poly);
    }
    else
    {
        nmod_poly_fit_length(poly, n);
        poly->length = n;
        _nmod_poly_interpolate_nmod_vec_barycentric(poly->coeffs,
            xs, ys, n, poly->mod);
        _nmod_poly_normalise(poly);
    }
}
