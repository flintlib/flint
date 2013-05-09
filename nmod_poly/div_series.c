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
_nmod_poly_div_series(mp_ptr Q, mp_srcptr A, mp_srcptr B, 
                                             len_t n, nmod_t mod)
{
    mp_ptr Binv = _nmod_vec_init(n);

    _nmod_poly_inv_series(Binv, B, n, mod);
    _nmod_poly_mullow(Q, Binv, n, A, n, n, mod);

    _nmod_vec_clear(Binv);
}

void
nmod_poly_div_series(nmod_poly_t Q, const nmod_poly_t A, 
                                    const nmod_poly_t B, len_t n)
{
    mp_ptr A_coeffs, B_coeffs, Q_coeffs;
    nmod_poly_t t1;
    len_t Alen, Blen;
    
    Blen = B->length;

    if (n == 0 || Blen == 0 || B->coeffs[0] == 0)
    {
        printf("Exception (nmod_poly_div_series). Division by zero.\n");
        abort();
    }
    
    Alen = A->length;

    if (Alen < n)
    {
        A_coeffs = _nmod_vec_init(n);
        flint_mpn_copyi(A_coeffs, A->coeffs, Alen);
        flint_mpn_zero(A_coeffs + Alen, n - Alen);
    }
    else
        A_coeffs = A->coeffs;

    if (Blen < n)
    {
        B_coeffs = _nmod_vec_init(n);
        flint_mpn_copyi(B_coeffs, B->coeffs, Blen);
        flint_mpn_zero(B_coeffs + Blen, n - Blen);
    }
    else
        B_coeffs = B->coeffs;

    if ((Q == A || Q == B) && Q->length >= n)
    {
        nmod_poly_init2(t1, Q->mod.n, n);
        Q_coeffs = t1->coeffs;
    }
    else
    {
        nmod_poly_fit_length(Q, n);
        Q_coeffs = Q->coeffs;
    }

    _nmod_poly_div_series(Q_coeffs, A_coeffs, B_coeffs, n, Q->mod);

    if ((Q == A || Q == B) && Q->length >= n)
    {
        nmod_poly_swap(Q, t1);
        nmod_poly_clear(t1);
    }
    
    Q->length = n;

    if (Alen < n)
        _nmod_vec_clear(A_coeffs);

    if (Blen < n)
        _nmod_vec_clear(B_coeffs);

    _nmod_poly_normalise(Q);
}
