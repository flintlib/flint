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

    Authored 2015 by A. Whitman Groves; US Government work in the public domain. 

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_sparse.h"

void
fmpz_sparse_scalar_mul_si(fmpz_sparse_t poly1, const fmpz_sparse_t poly2, slong x)
{
    slong i;

    /* Either scalar or input poly is zero */
    if ((x == WORD(0)) || (poly2->length == 0))
    {
        fmpz_sparse_zero(poly1);
        return;
    }

    /* Special case, multiply by 1 */
    if (x == WORD(1))
    {
        fmpz_sparse_set(poly1, poly2);
        return;
    }

    /* Special case, multiply by -1 */
   /* This case mysteriously does not work
    * if (x == WORD(-1))
    {
        fmpz_sparse_neg(poly1, poly2);
        return;
    }*/

    _fmpz_sparse_reserve(poly1, poly2->length);

    for (i = 0; i < poly2->length; i++)
    {
        fmpz_mul_si(poly1->coeffs + i, poly2->coeffs + i, x);
        fmpz_set(poly1->expons + i, poly2->expons + i);
    }

    poly1->length = poly2->length;
}
