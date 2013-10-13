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
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fq_poly.h"

void _fq_poly_div_newton_preinv (fq_struct *Q, const fq_struct *A, slong lenA,
                                 const fq_struct* B, slong lenB,
                                 const fq_struct* Binv, slong lenBinv,
                                 const fq_ctx_t ctx)
{
    const slong lenQ = lenA - lenB + 1;
    fq_struct *Arev, *Brev;

    Arev = _fq_vec_init(2*lenQ, ctx);
    Brev= Arev + lenQ;

    _fq_poly_reverse(Arev, A + (lenA - lenQ), lenQ, lenQ, ctx);

    Brev = Binv;
    _fq_vec_zero(Brev + lenBinv, lenQ - lenBinv, ctx);

    _fq_poly_mullow(Q, Arev, lenQ, Binv, lenQ, lenQ, ctx);

    _fq_poly_reverse(Q, Q, lenQ, lenQ, ctx);

    _fq_vec_clear(Arev, 2*lenQ);
}

void fq_poly_div_newton_preinv (fq_poly_t Q, const fq_poly_t A,
                                const fq_poly_t B, const fq_poly_t Binv,
                                const fq_ctx_t ctx)
{
    const slong lenA = A->length,
        lenB = B->length,
        lenQ = lenA - lenB + 1,
        lenBinv= Binv->length;

    fq_struct* q;

    if (lenB == 0)
    {
        printf("Exception (fq_poly_div_newton). Division by zero.\n");
        abort();
    }

    if (lenA < lenB)
    {
        fq_poly_zero(Q, ctx);
        return;
    }

    if (Q == A || Q == B || Q == Binv)
    {
        q = _fq_vec_init(lenQ, ctx);
    }
    else
    {
        fq_poly_fit_length(Q, lenQ, ctx);
        q = Q->coeffs;
    }

    _fq_poly_div_newton_preinv(q, A->coeffs, lenA, B->coeffs, lenB,
                               Binv->coeffs, lenBinv, ctx);

    if (Q == A || Q == B || Q == Binv)
    {
        flint_free(Q->coeffs);
        Q->coeffs = q;
        Q->alloc  = lenQ;
    }
    Q->length = lenQ;
}
