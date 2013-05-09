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

#define FLINT_REVERSE_NEWTON_CUTOFF 15

void
_nmod_poly_revert_series_newton(mp_ptr Qinv, mp_srcptr Q, len_t n, nmod_t mod)
{
    len_t *a, i, k;
    mp_ptr T, U, V;

    if (n >= 1) Qinv[0] = 0UL;
    if (n >= 2) Qinv[1] = n_invmod(Q[1], mod.n);
    if (n <= 2)
        return;

    T = _nmod_vec_init(n);
    U = _nmod_vec_init(n);
    V = _nmod_vec_init(n);

    k = n;
    for (i = 1; (1L << i) < k; i++);
    a = (len_t *) flint_malloc(i * sizeof(len_t));
    a[i = 0] = k;
    while (k >= FLINT_REVERSE_NEWTON_CUTOFF)
        a[++i] = (k = (k + 1) / 2);

    _nmod_poly_revert_series_lagrange(Qinv, Q, k, mod);
    _nmod_vec_zero(Qinv + k, n - k);

    for (i--; i >= 0; i--)
    {
        k = a[i];
        _nmod_poly_compose_series(T, Q, k, Qinv, k, k, mod);
        _nmod_poly_derivative(U, T, k, mod); U[k - 1] = 0UL;
        T[1] = 0UL;
        _nmod_poly_div_series(V, T, U, k, mod);
        _nmod_poly_derivative(T, Qinv, k, mod);
        _nmod_poly_mullow(U, V, k, T, k, k, mod);
        _nmod_vec_sub(Qinv, Qinv, U, k, mod);
    }

    flint_free(a);
    _nmod_vec_clear(T);
    _nmod_vec_clear(U);
    _nmod_vec_clear(V);
}

void
nmod_poly_revert_series_newton(nmod_poly_t Qinv, 
                                 const nmod_poly_t Q, len_t n)
{
    mp_ptr Qinv_coeffs, Q_coeffs;
    nmod_poly_t t1;
    len_t Qlen;
    
    Qlen = Q->length;

    if (Qlen < 2 || Q->coeffs[0] != 0 || Q->coeffs[1] == 0)
    {
        printf("Exception (nmod_poly_revert_series_newton). Input must have \n"
               "zero constant and an invertible coefficient of x^1.\n");
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

    _nmod_poly_revert_series_newton(Qinv_coeffs, Q_coeffs, n, Q->mod);

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
