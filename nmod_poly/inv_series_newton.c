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

#define NMOD_NEWTON_INVERSE_CUTOFF 400

void
_nmod_poly_inv_series_newton(mp_ptr Qinv, 
                                  mp_srcptr Q, long n, nmod_t mod)
{
    long m;
    mp_ptr prod, prod2;

    if (n < NMOD_NEWTON_INVERSE_CUTOFF)
    {
        _nmod_poly_inv_series_basecase(Qinv, Q, n, mod);
        return;
    }

    m = (n + 1)/2;
    prod = _nmod_vec_init(n);
    prod2 = _nmod_vec_init(n);

    _nmod_poly_inv_series_newton(Qinv, Q, m, mod);
    _nmod_poly_mullow(prod, Q, n, Qinv, m, n, mod);
    _nmod_poly_mullow(prod2 + m, Qinv, m, prod + m, n - m, n - m, mod);
    _nmod_vec_neg(Qinv + m, prod2 + m, n - m, mod);

    _nmod_vec_free(prod2);
    _nmod_vec_free(prod);
}

void
nmod_poly_inv_series_newton(nmod_poly_t Qinv, 
                                 const nmod_poly_t Q, long n)
{
    mp_ptr Qinv_coeffs, Q_coeffs;
    nmod_poly_t t1;
    long Qlen;
    
    Qlen = Q->length;

    if (n == 0 || Q->length == 0 || Q->coeffs[0] == 0)
    {
        printf("Exception: division by zero in nmod_poly_inv_series_newton\n");
        abort();
    }

    if (Qlen < n)
    {
        Q_coeffs = _nmod_vec_init(n);
        mpn_copyi(Q_coeffs, Q->coeffs, Qlen);
        mpn_zero(Q_coeffs + Qlen, n - Qlen);
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

    _nmod_poly_inv_series_newton(Qinv_coeffs, Q_coeffs, n, Q->mod);

    if (Q == Qinv && Qlen >= n)
    {
        nmod_poly_swap(Qinv, t1);
        nmod_poly_clear(t1);
    }
    
    Qinv->length = n;

    if (Qlen < n)
        _nmod_vec_free(Q_coeffs);

    _nmod_poly_normalise(Qinv);
}
