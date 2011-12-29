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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void _nmod_poly_divrem_q0(mp_ptr Q, mp_ptr R, 
                          mp_srcptr A, mp_srcptr B, long lenA, nmod_t mod)
{
    const mp_limb_t invL = (B[lenA-1] == 1) ? 1 : n_invmod(B[lenA-1], mod.n);

    if (lenA == 1)
    {
        _nmod_vec_scalar_mul_nmod(Q, A, lenA, invL, mod);
    }
    else
    {
        Q[0] = n_mulmod2_preinv(A[lenA-1], invL, mod.n, mod.ninv);

        _nmod_vec_scalar_mul_nmod(R, B, lenA - 1, Q[0], mod);
        _nmod_vec_sub(R, A, R, lenA - 1, mod);
    }
}

void nmod_poly_divrem_q0(nmod_poly_t Q, nmod_poly_t R, 
                         const nmod_poly_t A, const nmod_poly_t B)
{
    const long lenA = A->length, lenB = B->length;
    mp_ptr q, r;

    if (lenA != lenB || lenB == 0)
    {
        printf("Exception (nmod_poly_divrem_q0).  Expected lenA = lenB > 0.\n");
        abort();
    }

    if (Q == A || Q == B)
    {
        q = _nmod_vec_init(1);
    }
    else
    {
        nmod_poly_fit_length(Q, 1);
        q = Q->coeffs;
    }
    if (R == A || R == B)
    {
        r = _nmod_vec_init(lenB - 1);
    }
    else
    {
        nmod_poly_fit_length(R, lenB - 1);
        r = R->coeffs;
    }

    _nmod_poly_divrem_q0(q, r, A->coeffs, B->coeffs, lenB, A->mod);

    if (Q == A || Q == B)
    {
        _nmod_vec_clear(Q->coeffs);
        Q->coeffs = q;
        Q->alloc  = 1;
    }
    if (R == A || R == B)
    {
        _nmod_vec_clear(R->coeffs);
        R->coeffs = r;
        R->alloc  = lenB - 1;
    }
    Q->length = 1;
    R->length = lenB - 1;

    _nmod_poly_normalise(R);
}

