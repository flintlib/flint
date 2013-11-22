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

void _nmod_poly_div_newton_n_preinv (mp_ptr Q, mp_srcptr A, slong lenA,
                                     mp_srcptr B, slong lenB, mp_srcptr Binv,
                                     slong lenBinv, nmod_t mod)
{
    const slong lenQ = lenA - lenB + 1;
    mp_ptr Arev;

    Arev = _nmod_vec_init(lenQ);
    _nmod_poly_reverse(Arev, A + (lenA - lenQ), lenQ, lenQ);

    _nmod_poly_mullow(Q, Arev, lenQ, Binv, FLINT_MIN (lenQ,lenBinv), lenQ, mod);

    _nmod_poly_reverse(Q, Q, lenQ, lenQ);

    _nmod_vec_clear(Arev);
}

void nmod_poly_div_newton_n_preinv (nmod_poly_t Q, const nmod_poly_t A,
                                    const nmod_poly_t B, const nmod_poly_t Binv)
{
    const slong lenA = A->length, lenB = B->length, lenQ = lenA - lenB + 1,
               lenBinv = Binv->length;

    mp_ptr q;

    if (lenB == 0)
    {
        flint_printf("Exception (nmod_poly_div_newton_n_preinv). Division by zero.\n");
        abort();
    }

    if (lenA < lenB)
    {
        nmod_poly_zero(Q);
        return;
    }

    if (lenA > 2 * lenB - 2)
    {
        flint_printf ("Exception (nmod_poly_div_newton_n_preinv).\n");
    }

    if (Q == A || Q == B || Q == Binv)
    {
        q = flint_malloc(lenQ * sizeof(mp_limb_t));
    }
    else
    {
        nmod_poly_fit_length(Q, lenQ);
        q = Q->coeffs;
    }

    _nmod_poly_div_newton_n_preinv (q, A->coeffs, lenA, B->coeffs, lenB,
                                    Binv->coeffs, lenBinv, B->mod);

    if (Q == A || Q == B || Q == Binv)
    {
        flint_free(Q->coeffs);
        Q->coeffs = q;
        Q->alloc  = lenQ;
    }
    Q->length = lenQ;
}
