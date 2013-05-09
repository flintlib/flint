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

#include "nmod_poly.h"
#include "ulong_extras.h"

void
_nmod_poly_product_roots_nmod_vec(mp_ptr poly, mp_srcptr xs, len_t n, nmod_t mod)
{
    if (n == 0)
    {
        poly[0] = 1UL;
    }
    else if (n < 20)
    {
        len_t i, j;

        poly[n] = 1UL;
        poly[n - 1] = nmod_neg(xs[0], mod);

        for (i = 1; i < n; i++)
        {
            poly[n-i-1] = nmod_neg(n_mulmod2_preinv(poly[n-i], xs[i],
                mod.n, mod.ninv), mod);

            for (j = 0; j < i - 1; j++)
            {
                poly[n-i+j] = nmod_sub(poly[n-i+j],
                    n_mulmod2_preinv(poly[n-i+j+1], xs[i], mod.n, mod.ninv),
                        mod);
            }

            poly[n-1] = nmod_sub(poly[n-1], xs[i], mod);
        }
    }
    else
    {
        const len_t m = (n + 1) / 2;
        mp_ptr tmp;

        tmp = _nmod_vec_init(n + 2);

        _nmod_poly_product_roots_nmod_vec(tmp, xs, m, mod);
        _nmod_poly_product_roots_nmod_vec(tmp + m + 1, xs + m, n - m, mod);
        _nmod_poly_mul(poly, tmp, m + 1, tmp + m + 1, n - m + 1, mod);

        _nmod_vec_clear(tmp);
    }
}

void
nmod_poly_product_roots_nmod_vec(nmod_poly_t poly, mp_srcptr xs, len_t n)
{
    nmod_poly_fit_length(poly, n + 1);
    _nmod_poly_product_roots_nmod_vec(poly->coeffs, xs, n, poly->mod);
    poly->length = n + 1;
}
