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
#include "fmpz_mat.h"

void
fmpz_poly_mat_mul_KS(fmpz_poly_mat_t C, const fmpz_poly_mat_t A,
    const fmpz_poly_mat_t B)
{
    long i, j;
    long A_len, B_len;
    int signs;
    long A_bits, B_bits, bit_size;

    fmpz_mat_t AA, BB, CC;

    if (B->r == 0)
    {
        fmpz_poly_mat_zero(C);
        return;
    }

    A_len = fmpz_poly_mat_max_length(A);
    B_len = fmpz_poly_mat_max_length(B);

    A_bits = fmpz_poly_mat_max_bits(A);
    B_bits = fmpz_poly_mat_max_bits(B);

    signs = (A_bits < 0 || B_bits < 0);

    bit_size = FLINT_ABS(A_bits) + FLINT_ABS(B_bits) + signs;
    bit_size += FLINT_BIT_COUNT(FLINT_MIN(A_len, B_len));
    bit_size += FLINT_BIT_COUNT(B->r);

/*
    printf("A: BITS %ld LEN %ld\n", A_bits, A_len);
    printf("B: BITS %ld LEN %ld\n", B_bits, B_len);
    printf("bit_size: %ld\n", bit_size);
*/

    fmpz_mat_init(AA, A->r, A->c);
    fmpz_mat_init(BB, B->r, B->c);
    fmpz_mat_init(CC, C->r, C->c);

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_poly_bit_pack(fmpz_mat_entry(AA, i, j),
                               fmpz_poly_mat_entry(A, i, j), bit_size);

    for (i = 0; i < B->r; i++)
        for (j = 0; j < B->c; j++)
            fmpz_poly_bit_pack(fmpz_mat_entry(BB, i, j),
                               fmpz_poly_mat_entry(B, i, j), bit_size);

    fmpz_mat_mul(CC, AA, BB);

    for (i = 0; i < C->r; i++)
        for (j = 0; j < C->c; j++)
            if (signs)
                fmpz_poly_bit_unpack(fmpz_poly_mat_entry(C, i, j),
                    fmpz_mat_entry(CC, i, j), bit_size);
            else
                fmpz_poly_bit_unpack_unsigned(fmpz_poly_mat_entry(C, i, j),
                    fmpz_mat_entry(CC, i, j), bit_size);

    fmpz_mat_clear(AA);
    fmpz_mat_clear(BB);
    fmpz_mat_clear(CC);
}
