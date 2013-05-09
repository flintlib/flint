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
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void _nmod_poly_div_newton(mp_ptr Q, mp_srcptr A, len_t lenA, 
                                     mp_srcptr B, len_t lenB, nmod_t mod)
{
    const len_t lenQ = lenA - lenB + 1;
    mp_ptr Arev, Brev;

    Arev = _nmod_vec_init(2 * lenQ);
    Brev = Arev + lenQ;

    _nmod_poly_reverse(Arev, A + (lenA - lenQ), lenQ, lenQ);

    if (lenB >= lenQ)
    {
        _nmod_poly_reverse(Brev, B + (lenB - lenQ), lenQ, lenQ);
    }
    else
    {
        _nmod_poly_reverse(Brev, B, lenB, lenB);
        flint_mpn_zero(Brev + lenB, lenQ - lenB);
    }

    _nmod_poly_div_series(Q, Arev, Brev, lenQ, mod);

    _nmod_poly_reverse(Q, Q, lenQ, lenQ);

    _nmod_vec_clear(Arev);
}

void nmod_poly_div_newton(nmod_poly_t Q, const nmod_poly_t A,
                                         const nmod_poly_t B)
{
    const len_t lenA = A->length, lenB = B->length, lenQ = lenA - lenB + 1;

    mp_ptr q;

    if (lenB == 0)
    {
        printf("Exception (nmod_poly_div_newton). Division by zero.\n");
        abort();
    }

    if (lenA < lenB)
    {
        nmod_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        q = flint_malloc(lenQ * sizeof(mp_limb_t));
    }
    else
    {
        nmod_poly_fit_length(Q, lenQ);
        q = Q->coeffs;
    }

    _nmod_poly_div_newton(q, A->coeffs, lenA, B->coeffs, lenB, B->mod);

    if (Q == A || Q == B)
    {
        flint_free(Q->coeffs);
        Q->coeffs = q;
        Q->alloc  = lenQ;
    }
    Q->length = lenQ;
}
