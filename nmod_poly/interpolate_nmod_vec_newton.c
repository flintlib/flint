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


static void
_interpolate_newton(mp_ptr ys, mp_srcptr xs, len_t n, nmod_t mod)
{
    mp_limb_t p, q, t;
    len_t i, j;

    for (i = 1; i < n; i++)
    {
        t = ys[i - 1];

        for (j = i; j < n; j++)
        {
            p = nmod_sub(ys[j], t, mod);
            q = nmod_sub(xs[j], xs[j - i], mod);
            t = ys[j];
            q = n_invmod(q, mod.n);
            ys[j] = n_mulmod2_preinv(p, q, mod.n, mod.ninv);
        }
    }
}

static void
_newton_to_monomial(mp_ptr ys, mp_srcptr xs, len_t n, nmod_t mod)
{
    mp_limb_t t;
    len_t i, j;

    for (i = n - 2; i >= 0; i--)
    {
        t = ys[i];
        ys[i] = ys[i + 1];

        for (j = i + 1; j < n - 1; j++)
        {
            ys[j] = nmod_sub(ys[j + 1],
                n_mulmod2_preinv(ys[j], xs[i], mod.n, mod.ninv), mod);
        }

        ys[n - 1] = nmod_sub(t,
            n_mulmod2_preinv(ys[n - 1], xs[i], mod.n, mod.ninv), mod);
    }

    _nmod_poly_reverse(ys, ys, n, n);
}

void
_nmod_poly_interpolate_nmod_vec_newton(mp_ptr poly, mp_srcptr xs,
    mp_srcptr ys, len_t n, nmod_t mod)
{
    if (n == 1)
    {
        poly[0] = ys[0];
    }
    else
    {
        _nmod_vec_set(poly, ys, n);
        _interpolate_newton(poly, xs, n, mod);
        while (n > 0 && !poly[n-1]) n--;
        _newton_to_monomial(poly, xs, n, mod);
    }
}

void
nmod_poly_interpolate_nmod_vec_newton(nmod_poly_t poly,
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
        _nmod_poly_interpolate_nmod_vec_newton(poly->coeffs,
            xs, ys, n, poly->mod);
        _nmod_poly_normalise(poly);
    }
}
