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
_fmpz_poly_mullow(fmpz * res, const fmpz * poly1, len_t len1, 
                                const fmpz * poly2, len_t len2, len_t n)
{
    mp_size_t limbs1, limbs2;

    if (len2 < 7 || n < 7)
    {
        _fmpz_poly_mullow_classical(res, poly1, len1, poly2, len2, n);
        return;
    }

    limbs1 = _fmpz_vec_max_limbs(poly1, len1);
    limbs2 = _fmpz_vec_max_limbs(poly2, len2);

    if (n < 16 && (limbs1 > 12 || limbs2 > 12))
    {
        int clear = 0, i;
        fmpz *copy1, *copy2;

        if (len1 >= n)
            copy1 = (fmpz *) poly1;
        else
        {
            copy1 = (fmpz *) flint_malloc(n * sizeof(fmpz));
            for (i = 0; i < len1; i++)
                copy1[i] = poly1[i];
            flint_mpn_zero((mp_ptr) copy1 + len1, n - len1);
            clear |= 1;
        }

        if (len2 >= n)
            copy2 = (fmpz *) poly2;
        else
        {
            copy2 = (fmpz *) flint_malloc(n * sizeof(fmpz));
            for (i = 0; i < len2; i++)
                copy2[i] = poly2[i];
            flint_mpn_zero((mp_ptr) copy2 + len2, n - len2);
            clear |= 2;
        }

        _fmpz_poly_mullow_karatsuba_n(res, copy1, copy2, n);

        if (clear & 1)
            flint_free(copy1);
        if (clear & 2)
            flint_free(copy2);
    }
    else if (limbs1 + limbs2 <= 8)
        _fmpz_poly_mullow_KS(res, poly1, len1, poly2, len2, n);
    else if ((limbs1+limbs2)/2048 > len1 + len2)
        _fmpz_poly_mullow_KS(res, poly1, len1, poly2, len2, n);
    else if ((limbs1 + limbs2)*FLINT_BITS*4 < len1 + len2)
        _fmpz_poly_mullow_KS(res, poly1, len1, poly2, len2, n);
    else
        _fmpz_poly_mullow_SS(res, poly1, len1, poly2, len2, n);
}

void
fmpz_poly_mullow(fmpz_poly_t res,
                   const fmpz_poly_t poly1, const fmpz_poly_t poly2,
                   len_t n)
{
    const len_t len1 = poly1->length;
    const len_t len2 = poly2->length;

    if (len1 == 0 || len2 == 0 || n == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        fmpz_poly_mullow(t, poly1, poly2, n);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
        return;
    }

    n = FLINT_MIN(n, len1 + len2 - 1);

    fmpz_poly_fit_length(res, n);
    if (len1 >= len2)
        _fmpz_poly_mullow(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, n);
    else
        _fmpz_poly_mullow(res->coeffs, poly2->coeffs, len2, poly1->coeffs, len1, n);
    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);
}
