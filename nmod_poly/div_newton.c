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
_nmod_poly_div_newton(mp_ptr Q, mp_srcptr A, long Alen, 
                      mp_srcptr B, long Blen, nmod_t mod)
{
    long len = Alen - Blen + 1;
    long i;
    mp_ptr Arev, Brev;

    Arev = nmod_vec_init(2*len);
    Brev = Arev + len;

    if (Alen >= len)
        _nmod_poly_reverse(Arev, A + Alen - len, len, len);
    else
    {
        _nmod_poly_reverse(Arev, A, Alen, Alen);
        mpn_zero(Arev + Alen, len - Alen);
    }

    if (Blen >= len)
        _nmod_poly_reverse(Brev, B + Blen - len, len, len);
    else
    {
        _nmod_poly_reverse(Brev, B, Blen, Blen);
        mpn_zero(Brev + Blen, len - Blen);
    }

    _nmod_poly_div_series(Q, Arev, Brev, len, mod);
    
    _nmod_poly_reverse(Q, Q, len, len);
    nmod_vec_free(Arev);
}

void
nmod_poly_div_newton(nmod_poly_t Q, const nmod_poly_t A,
                                            const nmod_poly_t B)
{
    mp_ptr Q_coeffs;
    nmod_poly_t t1;
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

    _nmod_poly_div_newton(Q_coeffs, A->coeffs, Alen,
                               B->coeffs, Blen, B->mod);

    if (Q == A || Q == B)
    {
        nmod_poly_swap(Q, t1);
        nmod_poly_clear(t1);
    }
    
    Q->length = Alen - Blen + 1;

    _nmod_poly_normalise(Q);
}
