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
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_mullow_n(fmpz * res, const fmpz * poly1, const fmpz * poly2, long n)
{
    mp_size_t limbs1, limbs2;

    if (n < 7)
    {
        _fmpz_poly_mullow_classical(res, poly1, n, poly2, n, n);
        return;
    }

    limbs1 = _fmpz_vec_max_limbs(poly1, n);
    limbs2 = _fmpz_vec_max_limbs(poly2, n);

    if (n < 25 && (limbs1 > 12 || limbs2 > 12))
        _fmpz_poly_mullow_karatsuba_n(res, poly1, poly2, n);
    else
        _fmpz_poly_mullow_KS(res, poly1, n, poly2, n, n);
}

void
fmpz_poly_mullow_n(fmpz_poly_t res,
                   const fmpz_poly_t poly1, const fmpz_poly_t poly2,
                   long n)
{
    const long len1 = poly1->length;
    const long len2 = poly2->length;
    fmpz *copy1, *copy2;

    if (len1 == 0 || len2 == 0 || n == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        fmpz_poly_mullow_n(t, poly1, poly2, n);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
        return;
    }

    copy1 = poly1->coeffs;
    if (len1 < n)
    {
        long i;
        copy1 = (fmpz *) malloc(n * sizeof(fmpz));
        for (i = 0; i < len1; i++)
            copy1[i] = poly1->coeffs[i];
        for ( ; i < n; i++)
            copy1[i] = 0;
    }
    copy2 = (poly1 == poly2) ? copy1 : poly2->coeffs;
    if (poly1 != poly2 && len2 < n)
    {
        long i;
        copy2 = (fmpz *) malloc(n * sizeof(fmpz));
        for (i = 0; i < len2; i++)
            copy2[i] = poly2->coeffs[i];
        for ( ; i < n; i++)
            copy2[i] = 0;
    }

    fmpz_poly_fit_length(res, n);
    if (len1 >= len2)
        _fmpz_poly_mullow_n(res->coeffs, copy1, copy2, n);
    else
        _fmpz_poly_mullow_n(res->coeffs, copy2, copy1, n);
    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);

    if (copy1 != poly1->coeffs)
        free(copy1);
    if (poly1 != poly2 && copy2 != poly2->coeffs)
        free(copy2);
}
