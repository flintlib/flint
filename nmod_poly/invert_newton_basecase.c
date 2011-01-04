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
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void
_nmod_poly_invert_newton_basecase(mp_ptr Qinv, 
                                  mp_srcptr Q, long n, nmod_t mod)
{
    mp_ptr X2n;

    X2n = nmod_vec_init(2*n - 1);
    mpn_zero(X2n, 2*n - 2);
    X2n[2*n - 2] = 1;

    _nmod_poly_div_divconquer(Qinv, X2n, 2*n - 1, Q, n, mod);

    nmod_vec_free(X2n);
}

void
nmod_poly_invert_newton_basecase(nmod_poly_t Qinv, 
                                 const nmod_poly_t Q, long n)
{
    mp_ptr Qinv_coeffs, Q_coeffs;
    nmod_poly_t t1;
    long Qlen;
    
    Qlen = Q->length;

    if (Qlen < n || n == 0)
    {
        printf("Exception: division by zero in nmod_poly_invert_newton_basecase\n");
        abort();
    }

    Q_coeffs = Q->coeffs + Qlen - n;

    if (Q == Qinv)
    {
        nmod_poly_init2(t1, Q->mod.n, n);
        Qinv_coeffs = t1->coeffs;
    }
    else
    {
        nmod_poly_fit_length(Qinv, n);
        Qinv_coeffs = Qinv->coeffs;
    }

    _nmod_poly_invert_newton_basecase(Qinv_coeffs, Q_coeffs, n, Q->mod);

    if (Q == Qinv && Qlen >= n)
    {
        nmod_poly_swap(Qinv, t1);
        nmod_poly_clear(t1);
    }
    
    Qinv->length = n;

    _nmod_poly_normalise(Qinv);
}
