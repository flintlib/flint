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
    Copyright (C) 2013 Martin Lee

******************************************************************************/

#undef ulong
#define ulong ulongxx/* interferes with system includes */

#include <stdlib.h>

#undef ulong

#include <gmp.h>

#define ulong mp_limb_t

#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void _nmod_poly_divrem_newton_n_preinv (mp_ptr Q, mp_ptr R, mp_srcptr A,
                                        slong lenA, mp_srcptr B, slong lenB,
                                       mp_srcptr Binv, slong lenBinv, nmod_t mod)
{
    const slong lenQ = lenA - lenB + 1;

    if (lenA == lenB + 1)
    {
        _nmod_poly_divrem_q1 (Q, R, A, lenA, B, lenB, mod);
        return;
    }

    _nmod_poly_div_newton_n_preinv (Q, A, lenA, B, lenB, Binv, lenBinv, mod);

    if (lenB > 1)
    {
        if (lenQ >= lenB - 1)
            _nmod_poly_mullow(R, Q, lenQ, B, lenB - 1, lenB - 1, mod);
        else
            _nmod_poly_mullow(R, B, lenB - 1, Q, lenQ, lenB - 1, mod);

        _nmod_vec_sub(R, A, R, lenB - 1, mod);
    }
}

void nmod_poly_divrem_newton_n_preinv(nmod_poly_t Q, nmod_poly_t R,
                                      const nmod_poly_t A, const nmod_poly_t B,
                                      const nmod_poly_t Binv)
{
    const slong lenA = A->length, lenB = B->length, lenBinv = Binv->length;
    mp_ptr q, r;

    if (lenB == 0)
    {
        flint_printf("Exception (nmod_poly_divrem_newton_n_preinv). Division by zero.\n");
        abort();
    }

    if (lenA < lenB)
    {
        nmod_poly_set(R, A);
        nmod_poly_zero(Q);
        return;
    }

    if (lenA > 2 * lenB - 2)
    {
        flint_printf ("Exception (nmod_poly_divrem_newton_n_preinv).\n");
    }

    if (Q == A || Q == B || Q == Binv)
    {
        q = _nmod_vec_init(lenA - lenB + 1);
    }
    else
    {
        nmod_poly_fit_length(Q, lenA - lenB + 1);
        q = Q->coeffs;
    }
    if (R == A || R == B || R == Binv)
    {
        r = _nmod_vec_init(lenB - 1);
    }
    else
    {
        nmod_poly_fit_length(R, lenB - 1);
        r = R->coeffs;
    }

    _nmod_poly_divrem_newton_n_preinv (q, r, A->coeffs, lenA,
                                       B->coeffs, lenB, Binv->coeffs,
                                       lenBinv, B->mod);

    if (Q == A || Q == B || Q == Binv)
    {
        _nmod_vec_clear(Q->coeffs);
        Q->coeffs = q;
        Q->alloc  = lenA - lenB + 1;
    }
    if (R == A || R == B || R == Binv)
    {
        _nmod_vec_clear(R->coeffs);
        R->coeffs = r;
        R->alloc  = lenB - 1;
    }
    Q->length = lenA - lenB + 1;
    R->length = lenB - 1;

    _nmod_poly_normalise(R);
}
