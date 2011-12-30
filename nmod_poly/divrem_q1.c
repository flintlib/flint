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

void _nmod_poly_divrem_q1(mp_ptr Q, mp_ptr R, 
                          mp_srcptr A, long lenA, mp_srcptr B, long lenB,
                          nmod_t mod)
{
    const mp_limb_t invL = (B[lenB-1] == 1) ? 1 : n_invmod(B[lenB-1], mod.n);

    if (lenB == 1)
    {
        _nmod_vec_scalar_mul_nmod(Q, A, lenA, invL, mod);
    }
    else
    {
        mp_limb_t t;

        Q[1] = n_mulmod2_preinv(A[lenA-1], invL, mod.n, mod.ninv);
        t = n_mulmod2_preinv(Q[1], B[lenB-2], mod.n, mod.ninv);
        t = n_submod(A[lenA-2], t, mod.n);
        Q[0] = n_mulmod2_preinv(t, invL, mod.n, mod.ninv);

        if (FLINT_BITS + 2 <= 2 * mod.norm)
        {
            mpn_mul_1(R, B, lenB - 1, Q[0]);
            if (lenB > 2) 
                mpn_addmul_1(R + 1, B, lenB - 2, Q[1]);
            _nmod_vec_reduce(R, R, lenB - 1, mod);
        }
        else
        {
            _nmod_vec_scalar_mul_nmod(R, B, lenB - 1, Q[0], mod);
            if (lenB > 2)
                _nmod_vec_scalar_addmul_nmod(R + 1, B, lenB - 2, Q[1], mod);
        }

        _nmod_vec_sub(R, A, R, lenB - 1, mod);
    }
}

void nmod_poly_divrem_q1(nmod_poly_t Q, nmod_poly_t R, 
                         const nmod_poly_t A, const nmod_poly_t B)
{
    const long lenA = A->length, lenB = B->length;
    mp_ptr q, r;

    if (lenB == 0)
    {
        printf("Exception: division by zero in nmod_poly_divrem_q1\n");
        abort();
    }
    if (lenA < lenB)
    {
        nmod_poly_set(R, A);
        nmod_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        q = _nmod_vec_init(lenA - lenB + 1);
    }
    else
    {
        nmod_poly_fit_length(Q, lenA - lenB + 1);
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

    _nmod_poly_divrem_q1(q, r, A->coeffs, lenA, B->coeffs, lenB, A->mod);

    if (Q == A || Q == B)
    {
        _nmod_vec_clear(Q->coeffs);
        Q->coeffs = q;
        Q->alloc  = lenA - lenB + 1;
    }
    if (R == A || R == B)
    {
        _nmod_vec_clear(R->coeffs);
        R->coeffs = r;
        R->alloc  = lenB - 1;
    }
    Q->length = lenA - lenB + 1;
    R->length = lenB - 1;

    _nmod_poly_normalise(R);
}

