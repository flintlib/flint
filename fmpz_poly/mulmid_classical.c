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

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

/* Assumes poly1 and poly2 are not length 0 and len1 >= len2. */
void
_fmpz_poly_mulmid_classical(fmpz * res, const fmpz * poly1,
                            long len1, const fmpz * poly2, long len2)
{
    if ((len1 == 1) && (len2 == 1)) /* Special case if the length of both inputs is 1 */
    {
        fmpz_mul(res, poly1, poly2);
    }
    else                        /* Ordinary case */
    {
        long i;

        /* Set res[i] = poly1[i]*poly2[0]  */
        _fmpz_vec_scalar_mul_fmpz(res, poly1 + len2 - 1, len1 - len2 + 1,
                                  poly2);

        /* out[i+j] += in1[i]*in2[j] */
        for (i = 0; i < len2 - 1; i++)
            _fmpz_vec_scalar_addmul_fmpz(res, poly2 + len2 - i - 1,
                                         FLINT_MIN(i + 1, len1 - len2 + 1),
                                         poly1 + i);
        for (; i < len1 - 1; i++)
            _fmpz_vec_scalar_addmul_fmpz(res + i - len2 + 2, poly2 + 1,
                                         FLINT_MIN(len2 - 1, len1 - i - 1),
                                         poly1 + i);
    }
}


void
fmpz_poly_mulmid_classical(fmpz_poly_t res,
                           const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    long len_out;

    if (poly1->length == 0 || poly2->length == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    len_out = poly1->length - poly2->length + 1;

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t temp;
        fmpz_poly_init2(temp, len_out);
        _fmpz_poly_mulmid_classical(temp->coeffs, poly1->coeffs, poly1->length,
                                    poly2->coeffs, poly2->length);
        fmpz_poly_swap(res, temp);
        fmpz_poly_clear(temp);
    }
    else
    {
        fmpz_poly_fit_length(res, len_out);
        _fmpz_poly_mulmid_classical(res->coeffs, poly1->coeffs, poly1->length,
                                    poly2->coeffs, poly2->length);
    }

    _fmpz_poly_set_length(res, len_out);
    _fmpz_poly_normalise(res);
}
