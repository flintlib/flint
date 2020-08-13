/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"
#include "fmpz_mat.h"

void
fmpz_poly_mat_sqr_KS(fmpz_poly_mat_t B, const fmpz_poly_mat_t A)
{
    fmpz_mat_t AA, BB;
    slong i, j, n;
    slong A_len;
    int signs;
    slong A_bits, bit_size;

    n = A->r;

    if (n == 0)
    {
        fmpz_poly_mat_zero(B);
        return;
    }

    A_len = fmpz_poly_mat_max_length(A);
    A_bits = fmpz_poly_mat_max_bits(A);

    signs = A_bits < 0;

    bit_size = 2 * FLINT_ABS(A_bits) + signs;
    bit_size += FLINT_BIT_COUNT(A_len);
    bit_size += FLINT_BIT_COUNT(n);

    fmpz_mat_init(AA, n, n);
    fmpz_mat_init(BB, n, n);

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            fmpz_poly_bit_pack(fmpz_mat_entry(AA, i, j),
                               fmpz_poly_mat_entry(A, i, j), bit_size);

    /* Should use fmpz_mat_sqr */
    fmpz_mat_mul(BB, AA, AA);

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            if (signs)
                fmpz_poly_bit_unpack(fmpz_poly_mat_entry(B, i, j),
                    fmpz_mat_entry(BB, i, j), bit_size);
            else
                fmpz_poly_bit_unpack_unsigned(fmpz_poly_mat_entry(B, i, j),
                    fmpz_mat_entry(BB, i, j), bit_size);

    fmpz_mat_clear(AA);
    fmpz_mat_clear(BB);
}
