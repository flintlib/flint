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

    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"


void
_nmod_poly_sqrt_series(mp_ptr g, mp_srcptr h, long n, nmod_t mod)
{
    mp_ptr t = _nmod_vec_init(n);
    _nmod_poly_invsqrt_series(t, h, n, mod);
    _nmod_poly_mullow(g, t, n, h, n, n, mod);
    _nmod_vec_clear(t);
}

void
nmod_poly_sqrt_series(nmod_poly_t g, 
                                 const nmod_poly_t h, long n)
{
    mp_ptr g_coeffs, h_coeffs;
    nmod_poly_t t1;
    long hlen;
    
    hlen = h->length;

    if (n == 0)
    {
        printf("Exception (nmod_poly_sqrt_series). Division by zero.\n");
        abort();
    }

    if (h->length == 0 || h->coeffs[0] != 1UL)
    {
        printf("Exception (nmod_poly_sqrt_series). Requires constant term 1.\n");
        abort();
    }

    if (hlen < n)
    {
        h_coeffs = _nmod_vec_init(n);
        flint_mpn_copyi(h_coeffs, h->coeffs, hlen);
        flint_mpn_zero(h_coeffs + hlen, n - hlen);
    }
    else
        h_coeffs = h->coeffs;

    if (h == g && hlen >= n)
    {
        nmod_poly_init2(t1, h->mod.n, n);
        g_coeffs = t1->coeffs;
    }
    else
    {
        nmod_poly_fit_length(g, n);
        g_coeffs = g->coeffs;
    }

    _nmod_poly_sqrt_series(g_coeffs, h_coeffs, n, h->mod);

    if (h == g && hlen >= n)
    {
        nmod_poly_swap(g, t1);
        nmod_poly_clear(t1);
    }
    
    g->length = n;

    if (hlen < n)
        _nmod_vec_clear(h_coeffs);

    _nmod_poly_normalise(g);
}
