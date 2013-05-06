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
_nmod_poly_cosh_series(mp_ptr f, mp_srcptr h, long n, nmod_t mod)
{
    mp_ptr g = _nmod_vec_init(n);
    _nmod_poly_exp_expinv_series(f, g, h, n, mod);
    _nmod_vec_add(f, f, g, n, mod);
    _nmod_vec_scalar_mul_nmod(f, f, n, n_invmod(2UL, mod.n), mod);
    _nmod_vec_clear(g);
}

void
nmod_poly_cosh_series(nmod_poly_t g, const nmod_poly_t h, long n)
{
    mp_ptr g_coeffs, h_coeffs;
    nmod_poly_t t1;
    long h_len;
    
    h_len = h->length;

    if (h_len > 0 && h->coeffs[0] != 0UL)
    {
        printf("Exception (nmod_poly_cosh_series). Constant term != 0.\n");
        abort();
    }

    if (h_len == 1 || n < 2)
    {
        nmod_poly_zero(g);
        if (n > 0)
            nmod_poly_set_coeff_ui(g, 0, 1UL);
        return;
    }

    if (h_len < n)
    {
        h_coeffs = _nmod_vec_init(n);
        flint_mpn_copyi(h_coeffs, h->coeffs, h_len);
        flint_mpn_zero(h_coeffs + h_len, n - h_len);
    }
    else
        h_coeffs = h->coeffs;

    if (h == g && h_len >= n)
    {
        nmod_poly_init2(t1, h->mod.n, n);
        g_coeffs = t1->coeffs;
    }
    else
    {
        nmod_poly_fit_length(g, n);
        g_coeffs = g->coeffs;
    }

    _nmod_poly_cosh_series(g_coeffs, h_coeffs, n, h->mod);

    if (h == g && h_len >= n)
    {
        nmod_poly_swap(g, t1);
        nmod_poly_clear(t1);
    }
    
    g->length = n;

    if (h_len < n)
        _nmod_vec_clear(h_coeffs);

    _nmod_poly_normalise(g);
}
