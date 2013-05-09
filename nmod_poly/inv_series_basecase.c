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

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void
_nmod_poly_inv_series_basecase(mp_ptr Qinv, 
                                  mp_srcptr Q, len_t n, nmod_t mod)
{
    mp_ptr X2n, Qrev;

    X2n = _nmod_vec_init(2*n);
    Qrev = X2n + n;

    _nmod_poly_reverse(Qrev, Q, n, n);
    
    X2n[n - 1] = 1;
    flint_mpn_zero(X2n, n - 1);
    X2n -= (n - 1);

    _nmod_poly_div_divconquer(Qinv, X2n, 2*n - 1, Qrev, n, mod);

    _nmod_poly_reverse(Qinv, Qinv, n, n);
    
    _nmod_vec_clear(X2n + n - 1);
}

void
nmod_poly_inv_series_basecase(nmod_poly_t Qinv, 
                                 const nmod_poly_t Q, len_t n)
{
    mp_ptr Qinv_coeffs, Q_coeffs;
    nmod_poly_t t1;
    len_t Qlen;
    
    Qlen = Q->length;

    if (n == 0 || Q->length == 0 || Q->coeffs[0] == 0)
    {
        printf("Exception (nmod_poly_inv_series_basecase). Division by zero.\n");
        abort();
    }

    if (Qlen < n)
    {
        Q_coeffs = _nmod_vec_init(n);
        flint_mpn_copyi(Q_coeffs, Q->coeffs, Qlen);
        flint_mpn_zero(Q_coeffs + Qlen, n - Qlen);
    }
    else
        Q_coeffs = Q->coeffs;

    if (Q == Qinv && Qlen >= n)
    {
        nmod_poly_init2(t1, Q->mod.n, n);
        Qinv_coeffs = t1->coeffs;
    }
    else
    {
        nmod_poly_fit_length(Qinv, n);
        Qinv_coeffs = Qinv->coeffs;
    }

    _nmod_poly_inv_series_basecase(Qinv_coeffs, Q_coeffs, n, Q->mod);

    if (Q == Qinv && Qlen >= n)
    {
        nmod_poly_swap(Qinv, t1);
        nmod_poly_clear(t1);
    }
    
    Qinv->length = n;

    if (Qlen < n)
        _nmod_vec_clear(Q_coeffs);

    _nmod_poly_normalise(Qinv);
}
