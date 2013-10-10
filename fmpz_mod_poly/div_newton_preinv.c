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

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

void _fmpz_mod_poly_div_newton_preinv (fmpz *Q, const fmpz* A, slong lenA,
                                         const fmpz* B, slong lenB, const fmpz* Binv,
                                         slong lenBinv, const fmpz_t p)
{
    const slong lenQ = lenA - lenB + 1;
    fmpz *Arev, *Brev;

    Arev = _fmpz_vec_init(2*lenQ);
    Brev= Arev + lenQ;

    _fmpz_mod_poly_reverse(Arev, A + (lenA - lenQ), lenQ, lenQ);

    Brev = Binv;
    flint_mpn_zero ((mp_ptr) Brev + lenBinv, lenQ - lenBinv);

    _fmpz_mod_poly_mullow(Q, Arev, lenQ, Binv, lenQ, p, lenQ);

    _fmpz_mod_poly_reverse(Q, Q, lenQ, lenQ);

    _fmpz_vec_clear(Arev, 2*lenQ);
}

void fmpz_mod_poly_div_newton_preinv (fmpz_mod_poly_t Q, const fmpz_mod_poly_t A,
                                        const fmpz_mod_poly_t B, const fmpz_mod_poly_t Binv)
{
    const slong lenA = A->length,
        lenB = B->length,
        lenQ = lenA - lenB + 1,
        lenBinv= Binv->length;

    fmpz* q;

    if (lenB == 0)
    {
        printf("Exception (fmpz_mod_poly_div_newton). Division by zero.\n");
        abort();
    }

    if (lenA < lenB)
    {
        fmpz_mod_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B || Q == Binv)
    {
        q = _fmpz_vec_init(lenQ);
    }
    else
    {
        fmpz_mod_poly_fit_length(Q, lenQ);
        q = Q->coeffs;
    }

    _fmpz_mod_poly_div_newton_preinv(q, A->coeffs, lenA, B->coeffs, lenB,
                                       Binv->coeffs, lenBinv, &B->p);

    if (Q == A || Q == B || Q == Binv)
    {
        flint_free(Q->coeffs);
        Q->coeffs = q;
        Q->alloc  = lenQ;
    }
    Q->length = lenQ;
}
