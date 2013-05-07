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
    Copyright (C) 2010, 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void _fmpz_poly_sqrlow(fmpz * res, const fmpz * poly, long len, long n)
{
    mp_size_t limbs;

    if (n < 7)
    {
        _fmpz_poly_sqrlow_classical(res, poly, len, n);
        return;
    }

    limbs = _fmpz_vec_max_limbs(poly, len);

    if (n < 16 && limbs > 12)
    {
        int i;
        fmpz *copy;

        copy = flint_malloc(n * sizeof(fmpz));
        for (i = 0; i < len; i++)
            copy[i] = poly[i];
        flint_mpn_zero((mp_ptr) copy + len, n - len);

        _fmpz_poly_sqrlow_karatsuba_n(res, copy, n);

        flint_free(copy);
    }
    else if (limbs <= 4)
        _fmpz_poly_sqrlow_KS(res, poly, len, n);
    else if (limbs/2048 > len)
        _fmpz_poly_sqrlow_KS(res, poly, len, n);
    else if (limbs*FLINT_BITS*4 < len)
       _fmpz_poly_sqrlow_KS(res, poly, len, n);
    else
       _fmpz_poly_mullow_SS(res, poly, len, poly, len, n);
}

void fmpz_poly_sqrlow(fmpz_poly_t res, const fmpz_poly_t poly, long n)
{
    const long len = poly->length;

    if (len == 0 || n == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (res == poly)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        fmpz_poly_sqrlow(t, poly, n);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
        return;
    }

    n = FLINT_MIN(2 * len - 1, n);

    fmpz_poly_fit_length(res, n);
    _fmpz_poly_sqrlow(res->coeffs, poly->coeffs, len, n);
    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);
}
