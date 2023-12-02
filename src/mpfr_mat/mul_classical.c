/*
    Copyright (C) 2010,2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpfr_mat.h"

void
mpfr_mat_mul_classical(mpfr_mat_t C, const mpfr_mat_t A, const mpfr_mat_t B,
                       mpfr_rnd_t rnd)
{
    slong ar, bc, br;
    slong i, j, k;
    mpfr_t tmp;

    ar = A->r;
    br = B->r;
    bc = B->c;

    if (C == A || C == B)
    {
        mpfr_mat_t t;
        mpfr_mat_init(t, ar, bc, C->prec);
        mpfr_mat_mul_classical(t, A, B, rnd);
        mpfr_mat_swap_entrywise(C, t);
        mpfr_mat_clear(t);
        return;
    }

    if (C->r != ar || C->c != bc)
    {
        flint_throw(FLINT_ERROR, "(mpfr_mat_mul_classical): Incompatible dimensions.\n");
    }

    if (br == 0)
    {
        mpfr_mat_zero(C);
        return;
    }

    mpfr_init2(tmp, C->prec);

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            mpfr_mul(mpfr_mat_entry(C, i, j), mpfr_mat_entry(A, i, 0),
                     mpfr_mat_entry(B, 0, j), rnd);

            for (k = 1; k < br; k++)
            {
                mpfr_mul(tmp, mpfr_mat_entry(A, i, k), mpfr_mat_entry(B, k, j),
                         rnd);
                mpfr_add(mpfr_mat_entry(C, i, j), mpfr_mat_entry(C, i, j), tmp,
                         rnd);
            }
        }
    }

    mpfr_clear(tmp);
}
