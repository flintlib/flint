/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"


void
_nmod_poly_interpolate_nmod_vec_barycentric(mp_ptr poly,
                            mp_srcptr xs, mp_srcptr ys, len_t n, nmod_t mod)
{
    mp_ptr P, Q, w;
    len_t i, j;

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
        w[i] = 1UL;
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
                                    mp_srcptr xs, mp_srcptr ys, len_t n)
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
