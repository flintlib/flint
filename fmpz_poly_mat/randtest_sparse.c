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
#include "ulong_extras.h"

void
fmpz_poly_mat_randtest_sparse(fmpz_poly_mat_t A, flint_rand_t state, len_t len,
    mp_bitcnt_t bits, float density)
{
    len_t i, j;

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < A->c; j++)
        {
            if (n_randint(state, 1000) < density * 1000)
            {
                len_t l = n_randint(state, len + 1);
                l = FLINT_MAX(l, 1);
                fmpz_poly_randtest_not_zero(fmpz_poly_mat_entry(A, i, j),
                    state, l, bits);
            }
            else
            {
                fmpz_poly_zero(fmpz_poly_mat_entry(A, i, j));
            }
        }
    }
}
