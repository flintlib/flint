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
#include "nmod_poly.h"
#include "ulong_extras.h"

static void 
__nmod_poly_invsqrt_series_prealloc(mp_ptr g, 
                                    mp_srcptr h, mp_ptr t, mp_ptr u,
                                    len_t n, nmod_t mod)
{
    const int alloc = (t == NULL);
    const len_t m    = (n + 1) / 2;
    mp_limb_t c;

    if (n == 1)
    {
        g[0] = 1UL;
        return;
    }

    if (alloc)
    {
        t = _nmod_vec_init(n);
        u = _nmod_vec_init(n);
    }

    __nmod_poly_invsqrt_series_prealloc(g, h, t, u, m, mod);

    _nmod_vec_zero(g + m, n - m);

    _nmod_poly_mul(t, g, m, g, m, mod);
    if (2*m - 1 < n)
        t[n-1] = 0UL;

    _nmod_poly_mullow(u, t, n, g, n, n, mod);
    _nmod_poly_mullow(t, u, n, h, n, n, mod);

    c = n_invmod(mod.n - 2UL, mod.n);
    _nmod_vec_scalar_mul_nmod(g + m, t + m, n - m, c, mod);

    if (alloc)
    {
        _nmod_vec_clear(t);
        _nmod_vec_clear(u);
    }
}

void _nmod_poly_invsqrt_series(mp_ptr g, mp_srcptr h, len_t n, nmod_t mod)
{
    __nmod_poly_invsqrt_series_prealloc(g, h, NULL, NULL, n, mod);
}

void nmod_poly_invsqrt_series(nmod_poly_t g, const nmod_poly_t h, len_t n)
{
    const len_t hlen = h->length;
    mp_ptr g_coeffs, h_coeffs;
    nmod_poly_t t1;

    if (n == 0 || h->length == 0 || h->coeffs[0] == 0)
    {
        printf("Exception (nmod_poly_invsqrt). Division by zero.\n");
        abort();
    }

    if (h->coeffs[0] != 1UL)
    {
        printf("Exception (nmod_poly_invsqrt_series). Constant term != 1.\n");
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

    _nmod_poly_invsqrt_series(g_coeffs, h_coeffs, n, h->mod);

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
