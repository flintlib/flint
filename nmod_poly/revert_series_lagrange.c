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

   Copyright (C) 2010 Sebastian Pancratz
   Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void
_nmod_poly_revert_series_lagrange(mp_ptr Qinv, mp_srcptr Q, slong n, nmod_t mod)
{
    slong i;
    mp_ptr R, S, T, tmp;

    if (n >= 1) Qinv[0] = UWORD(0);
    if (n >= 2) Qinv[1] = n_invmod(Q[1], mod.n);
    if (n <= 2)
        return;

    R = _nmod_vec_init(n - 1);
    S = _nmod_vec_init(n - 1);
    T = _nmod_vec_init(n - 1);

    _nmod_poly_inv_series(R, Q + 1, n - 1, n - 1, mod);
    _nmod_vec_set(S, R, n - 1);

    for (i = 2; i < n; i++)
    {
        _nmod_poly_mullow(T, S, n - 1, R, n - 1, n - 1, mod);
        Qinv[i] = nmod_div(T[i - 1], i, mod);
        tmp = S; S = T; T = tmp;
    }

    _nmod_vec_clear(R);
    _nmod_vec_clear(S);
    _nmod_vec_clear(T);
}

void
nmod_poly_revert_series_lagrange(nmod_poly_t Qinv, 
                                 const nmod_poly_t Q, slong n)
{
    mp_ptr Qinv_coeffs, Q_coeffs;
    nmod_poly_t t1;
    slong Qlen;
    
    Qlen = Q->length;

    if (Qlen < 2 || Q->coeffs[0] != 0 || Q->coeffs[1] == 0)
    {
        flint_printf("Exception (nmod_poly_revert_series_lagrange). Input must \n"
               "have zero constant and an invertible coefficient of x^1.\n");
        flint_abort();
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

    _nmod_poly_revert_series_lagrange(Qinv_coeffs, Q_coeffs, n, Q->mod);

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
