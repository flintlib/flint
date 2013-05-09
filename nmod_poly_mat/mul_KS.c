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
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "fmpz.h"
#include "fmpz_mat.h"

void
nmod_poly_mat_mul_KS(nmod_poly_mat_t C, const nmod_poly_mat_t A,
    const nmod_poly_mat_t B)
{
    len_t i, j;
    len_t A_len, B_len;
    mp_bitcnt_t bit_size;

    fmpz_mat_t AA, BB, CC;

    if (B->r == 0)
    {
        nmod_poly_mat_zero(C);
        return;
    }

    A_len = nmod_poly_mat_max_length(A);
    B_len = nmod_poly_mat_max_length(B);

    bit_size = 2 * FLINT_BIT_COUNT(nmod_poly_mat_modulus(A));
    bit_size += FLINT_BIT_COUNT(FLINT_MIN(A_len, B_len));
    bit_size += FLINT_BIT_COUNT(B->r);

    fmpz_mat_init(AA, A->r, A->c);
    fmpz_mat_init(BB, B->r, B->c);
    fmpz_mat_init(CC, C->r, C->c);

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            nmod_poly_bit_pack(fmpz_mat_entry(AA, i, j),
                               nmod_poly_mat_entry(A, i, j), bit_size);

    for (i = 0; i < B->r; i++)
        for (j = 0; j < B->c; j++)
            nmod_poly_bit_pack(fmpz_mat_entry(BB, i, j),
                               nmod_poly_mat_entry(B, i, j), bit_size);

    fmpz_mat_mul(CC, AA, BB);

    for (i = 0; i < C->r; i++)
        for (j = 0; j < C->c; j++)
                nmod_poly_bit_unpack(nmod_poly_mat_entry(C, i, j),
                    fmpz_mat_entry(CC, i, j), bit_size);

    fmpz_mat_clear(AA);
    fmpz_mat_clear(BB);
    fmpz_mat_clear(CC);
}
