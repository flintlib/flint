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
_nmod_poly_divrem_newton(mp_ptr Q, mp_ptr R, mp_srcptr A, long Alen, 
                      mp_srcptr B, long Blen, nmod_t mod)
{
    long len = Alen - Blen + 1;
    
    _nmod_poly_div_newton(Q, A, Alen, B, Blen, mod);

    if (Blen > 1)
    {
        if (len >= Blen - 1)
            _nmod_poly_mullow_n(R, Q, len, B, Blen - 1, Blen - 1, mod);
        else
            _nmod_poly_mullow_n(R, B, Blen - 1, Q, len, Blen - 1, mod);

        _nmod_vec_sub(R, A, R, Blen - 1, mod);
    }
}

void
nmod_poly_divrem_newton(nmod_poly_t Q, nmod_poly_t R, 
                        const nmod_poly_t A, const nmod_poly_t B)
{
    mp_ptr Q_coeffs, R_coeffs;
    nmod_poly_t t1, t2;
    long Alen, Blen;

    Blen = B->length;

    if (Blen == 0)
    {
        printf("Exception: division by zero in nmod_poly_div_newton\n");
        abort();
    }

    Alen = A->length;

    if (Alen < Blen)
    {
        nmod_poly_set(R, A);
        nmod_poly_zero(Q);

        return;
    }

    if (Q == A || Q == B)
    {
        nmod_poly_init2_preinv(t1, B->mod.n, B->mod.ninv,
                               Alen - Blen + 1);
        Q_coeffs = t1->coeffs;
    }
    else
    {
        nmod_poly_fit_length(Q, Alen - Blen + 1);
        Q_coeffs = Q->coeffs;
    }

    if (R == A || R == B)
    {
        nmod_poly_init2_preinv(t2, B->mod.n, B->mod.ninv,
                                Blen - 1);
        R_coeffs = t2->coeffs;
    }
    else
    {
        nmod_poly_fit_length(R, Blen - 1);
        R_coeffs = R->coeffs;
    }

    _nmod_poly_divrem_newton(Q_coeffs, R_coeffs, A->coeffs, Alen,
                               B->coeffs, Blen, B->mod);

    if (Q == A || Q == B)
    {
        nmod_poly_swap(Q, t1);
        nmod_poly_clear(t1);
    }
    
    if (R == A || R == B)
    {
        nmod_poly_swap(R, t2);
        nmod_poly_clear(t2);
    }
    
    Q->length = Alen - Blen + 1;
    R->length = Blen - 1;

    _nmod_poly_normalise(Q);
    _nmod_poly_normalise(R);
}
