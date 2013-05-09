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

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

int _fmpz_poly_divides(fmpz * q, const fmpz * a, 
                       len_t len1, const fmpz * b, len_t len2)
{
    fmpz * r = _fmpz_vec_init(len1);

    _fmpz_poly_divrem(q, r, a, len1, b, len2);

    FMPZ_VEC_NORM(r, len1);

    _fmpz_vec_clear(r, len1);

    return (len1 == 0);
}

int fmpz_poly_divides(fmpz_poly_t q, const fmpz_poly_t a, const fmpz_poly_t b)
{
    if (fmpz_poly_is_zero(b))
    {
        printf("Exception (fmpz_poly_divides). Division by zero.\n");
        abort();
    }
    if (fmpz_poly_is_zero(a))
    {
        fmpz_poly_zero(q);
        return 1;
    }
    if (a->length < b->length)
    {
        return 0;
    }

    {
        const len_t lenQ = a->length - b->length + 1;
        int res;

        if (q == a || q == b)
        {
            fmpz_poly_t t;

            fmpz_poly_init2(t, lenQ);
            res = _fmpz_poly_divides(t->coeffs, a->coeffs, a->length, b->coeffs, b->length);
            _fmpz_poly_set_length(t, lenQ);
            _fmpz_poly_normalise(t);
            fmpz_poly_swap(q, t);
            fmpz_poly_clear(t);
        }
        else
        {
            fmpz_poly_fit_length(q, lenQ);
            res = _fmpz_poly_divides(q->coeffs, a->coeffs, a->length, b->coeffs, b->length);
            _fmpz_poly_set_length(q, lenQ);
            _fmpz_poly_normalise(q);
        }
        return res;
    }
}
