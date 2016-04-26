/*
    Copyright (C) 2010,2011 Fredrik Johansson
    Copyright (C) 2016 Aaditya Thakkar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void
fmpz_mat_mul(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong dim, m, n, k;

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
        slong ab, bb, bits;

        ab = fmpz_mat_max_bits(A);
        bb = fmpz_mat_max_bits(B);

        ab = FLINT_ABS(ab);
        bb = FLINT_ABS(bb);

        bits = ab + bb + FLINT_BIT_COUNT(n) + 1;

        if (5*(ab + bb) > dim * dim || (bits > FLINT_BITS - 3 && dim < 60))
        {
            if ((ab + bb) * dim < 17000)
            {
                fmpz_mat_mul_classical_inline(C, A, B);
            }
            else
            {
                if (dim > 75 && (ab + bb) > 650)
                {
                    _fmpz_mat_mul_multi_mod(C, A, B, bits);
                }
                else
                {
                    fmpz_mat_mul_strassen(C, A, B);
                }
            }
        }
        else
        {
            _fmpz_mat_mul_multi_mod(C, A, B, bits);
        }
    }
}
