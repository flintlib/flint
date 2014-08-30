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
   
******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

int
fmpz_poly_equal_trunc(const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n)
{
    slong i, len1, len2;

    if (poly1 == poly2)
        return 1;               /* same polynomial */

    if (n < 0)
       n = 0;

    len1 = FLINT_MIN(poly1->length, n);
    len2 = FLINT_MIN(poly2->length, n);

    if (len1 < len2)
    {
       for (i = len1; i < len2; i++)
       {
          if (!fmpz_is_zero(poly2->coeffs + i))
             return 0;
       }
    } else if (len2 < len1)
    {
       for (i = len2; i < len1; i++)
       {
          if (!fmpz_is_zero(poly1->coeffs + i))
             return 0;
       }
    }

    for (i = 0; i < FLINT_MIN(len1, len2); i++) /* check if coefficients the same */
        if (!fmpz_equal(poly1->coeffs + i, poly2->coeffs + i))
            return 0;

    return 1;
}
