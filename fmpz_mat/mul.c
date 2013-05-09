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

    Copyright (C) 2010,2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"


void
fmpz_mat_mul(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    len_t dim, m, n, k;

    m = A->r;
    n = A->c;
    k = B->c;

    if (C == A || C == B)
    {
        fmpz_mat_t t;
        fmpz_mat_init(t, m, k);
        fmpz_mat_mul(t, A, B);
        fmpz_mat_swap(C, t);
        fmpz_mat_clear(t);
        return;
    }

    dim = FLINT_MIN(FLINT_MIN(m, n), k);

    if (dim < 12)
    {
        /* The inline version only benefits from large n */
        if (n <= 2)
            fmpz_mat_mul_classical(C, A, B);
        else
            fmpz_mat_mul_classical_inline(C, A, B);
    }
    else
    {
        len_t ab, bb, bits;

        ab = fmpz_mat_max_bits(A);
        bb = fmpz_mat_max_bits(B);

        ab = FLINT_ABS(ab);
        bb = FLINT_ABS(bb);

        bits = ab + bb + FLINT_BIT_COUNT(n) + 1;

        if (5*(ab + bb) > dim * dim || (bits > FLINT_BITS - 3 && dim < 60))
        {
            fmpz_mat_mul_classical_inline(C, A, B);
        }
        else
        {
            _fmpz_mat_mul_multi_mod(C, A, B, bits);
        }
    }
}
