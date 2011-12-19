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
    Copyright (C) 2011 Sebastian Pancratz

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
    long Qlen = Alen - Blen + 1;
    mp_ptr Arev, Brev;

    Arev = _nmod_vec_init(2 * Qlen);
    Brev = Arev + Qlen;

    _nmod_poly_reverse(Arev, A + (Alen - Qlen), Qlen, Qlen);

    if (Blen >= Qlen)
    {
        _nmod_poly_reverse(Brev, B + (Blen - Qlen), Qlen, Qlen);
    }
    else
    {
        _nmod_poly_reverse(Brev, B, Blen, Blen);
        mpn_zero(Brev + Blen, Qlen - Blen);
    }

    _nmod_poly_div_series(Q, Arev, Brev, Qlen, mod);

    _nmod_poly_reverse(Q, Q, Qlen, Qlen);

    _nmod_vec_clear(Arev);
}

void
nmod_poly_div_newton(nmod_poly_t Q, const nmod_poly_t A,
                                    const nmod_poly_t B)
{
    const long Alen = A->length;
    const long Blen = B->length;
    const long Qlen = Alen - Blen + 1;

    mp_ptr q;

    if (Blen == 0)
    {
        printf("Exception: division by zero in nmod_poly_div_newton\n");
        abort();
    }

    if (Alen < Blen)
    {
        nmod_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        q = malloc(Qlen * sizeof(mp_limb_t));
    }
    else
    {
        nmod_poly_fit_length(Q, Qlen);
        q = Q->coeffs;
    }

    _nmod_poly_div_newton(q, A->coeffs, Alen, B->coeffs, Blen, B->mod);

    if (Q == A || Q == B)
    {
        free(Q->coeffs);
        Q->coeffs = q;
        Q->alloc  = Qlen;
        Q->length = Qlen;
    }
    else
    {
        Q->length = Qlen;
    }
}
