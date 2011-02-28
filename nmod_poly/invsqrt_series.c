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
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"


static void 
__nmod_poly_invsqrt_series_prealloc(mp_ptr g, 
                                      mp_srcptr h, mp_ptr t, mp_ptr u,
                                      long n, nmod_t mod)
{
    int alloc;
    long m;
    mp_limb_t c;

    if (n == 1)
    {
        g[0] = 1UL;
        return;
    }

    m = (n + 1)/2;

    alloc = (t == NULL);
    if (alloc)
    {
        t = _nmod_vec_init(n);
        u = _nmod_vec_init(n);
    }

    __nmod_poly_invsqrt_series_prealloc(g, h, t, u, m, mod);

    _nmod_vec_zero(g + m, n - m);
    _nmod_poly_mullow(t, g, n, g, n, n, mod);
    _nmod_poly_mullow(u, t, n, g, n, n, mod);
    _nmod_poly_mullow(t, u, n, h, n, n, mod);

    c = n_invmod(mod.n - 2UL, mod.n);
    _nmod_vec_scalar_mul_nmod(g + m, t + m, n - m, c, mod);

    if (alloc)
    {
        _nmod_vec_free(t);
        _nmod_vec_free(u);
    }
}


void
_nmod_poly_invsqrt_series(mp_ptr g, 
                                  mp_srcptr h, long n, nmod_t mod)
{
    __nmod_poly_invsqrt_series_prealloc(g, h, NULL, NULL, n, mod);
}

void
nmod_poly_invsqrt_series(nmod_poly_t g, 
                                 const nmod_poly_t h, long n)
{
    mp_ptr g_coeffs, h_coeffs;
    nmod_poly_t t1;
    long hlen;
    
    hlen = h->length;

    if (n == 0 || h->length == 0 || h->coeffs[0] == 0)
    {
        printf("Exception: division by zero in nmod_poly_invsqrt_series\n");
        abort();
    }

    if (h->coeffs[0] != 1UL)
    {
        printf("Exception: nmod_poly_invsqrt_series requires constant term 1\n");
        abort();
    }

    if (hlen < n)
    {
        h_coeffs = _nmod_vec_init(n);
        mpn_copyi(h_coeffs, h->coeffs, hlen);
        mpn_zero(h_coeffs + hlen, n - hlen);
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
        _nmod_vec_free(h_coeffs);

    _nmod_poly_normalise(g);
}
