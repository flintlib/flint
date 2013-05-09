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

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_rem_basecase(fmpz * R, const fmpz * A, len_t lenA,
                                  const fmpz * B, len_t lenB)
{
    const fmpz * leadB = B + (lenB - 1);
    fmpz_t q;

    fmpz_init(q);

    if (R != A)
        _fmpz_vec_set(R, A, lenA);

    for ( ; lenA >= lenB; lenA--)
    {
        if (fmpz_cmpabs(R + lenA - 1, leadB) >= 0)
        {
            fmpz_fdiv_q(q, R + lenA - 1, leadB);
            _fmpz_vec_scalar_submul_fmpz(R + lenA - lenB, B, lenB, q);
        }
    }

    fmpz_clear(q);
}

void
fmpz_poly_rem_basecase(fmpz_poly_t R,
                       const fmpz_poly_t A, const fmpz_poly_t B)
{
    len_t lenr;
    fmpz *r;
    
    if (B->length == 0)
    {
        printf("Exception (fmpz_poly_rem_basecase). Division by zero.\n");
        abort();
    }
    if (A->length < B->length)
    {
        fmpz_poly_set(R, A);
        return;
    }

    lenr = A->length;
    if (R == B)
        r = _fmpz_vec_init(lenr);
    else
    {
        fmpz_poly_fit_length(R, lenr);
        r = R->coeffs;
    }

    _fmpz_poly_rem_basecase(r, A->coeffs, A->length,
                               B->coeffs, B->length);

    if (R == B)
    {
        _fmpz_vec_clear(R->coeffs, R->alloc);
        R->coeffs = r;
        R->alloc = lenr;
        R->length = lenr;
    }
    else
        _fmpz_poly_set_length(R, lenr);

    _fmpz_poly_normalise(R);
}
