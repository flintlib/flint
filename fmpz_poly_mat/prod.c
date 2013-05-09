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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

static void
binary_splitting(fmpz_poly_mat_t P, fmpz_poly_mat_t * const factors,
                                                            len_t n1, len_t n2)
{
    if (n2 - n1 <= 0)
    {
        fmpz_poly_mat_one(P);
    }
    else if (n2 - n1 == 1)
    {
        fmpz_poly_mat_set(P, factors[n1]);
    }
    else if (n2 - n1 == 2)
    {
        fmpz_poly_mat_mul(P, factors[n1], factors[n1 + 1]);
    }
    else
    {
        fmpz_poly_mat_t P1, P2;
        len_t m = (n1 + n2) / 2;

        fmpz_poly_mat_init(P1, P->r, P->c);
        fmpz_poly_mat_init(P2, P->r, P->c);

        binary_splitting(P1, factors, n1, m);
        binary_splitting(P2, factors, m, n2);

        fmpz_poly_mat_mul(P, P1, P2);

        fmpz_poly_mat_clear(P1);
        fmpz_poly_mat_clear(P2);
    }
}

void
fmpz_poly_mat_prod(fmpz_poly_mat_t res,
                        fmpz_poly_mat_t * const factors, len_t n)
{
    binary_splitting(res, factors, 0, n);
}
