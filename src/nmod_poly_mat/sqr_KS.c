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
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "fmpz.h"
#include "fmpz_mat.h"

void
nmod_poly_mat_sqr_KS(nmod_poly_mat_t B, const nmod_poly_mat_t A)
{
    slong i, j, n;
    slong A_len;
    flint_bitcnt_t bit_size;
    fmpz_mat_t AA, BB;

    n = A->r;

    if (n == 0)
    {
        nmod_poly_mat_zero(B);
        return;
    }

    A_len = nmod_poly_mat_max_length(A);

    bit_size = 2 * FLINT_BIT_COUNT(nmod_poly_mat_modulus(A));
    bit_size += FLINT_BIT_COUNT(A_len);
    bit_size += FLINT_BIT_COUNT(n);

    fmpz_mat_init(AA, n, n);
    fmpz_mat_init(BB, n, n);

    for (i = 0; i < n; i++)
        for (j = 0; j < A->c; j++)
            nmod_poly_bit_pack(fmpz_mat_entry(AA, i, j),
                               nmod_poly_mat_entry(A, i, j), bit_size);

    /* Should use fmpz_mat_sqr */
    fmpz_mat_mul(BB, AA, AA);

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
                nmod_poly_bit_unpack(nmod_poly_mat_entry(B, i, j),
                    fmpz_mat_entry(BB, i, j), bit_size);

    fmpz_mat_clear(AA);
    fmpz_mat_clear(BB);
}
