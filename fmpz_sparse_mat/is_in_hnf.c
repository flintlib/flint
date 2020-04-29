/*
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_mat.h"

int fmpz_sparse_mat_is_in_hnf(const fmpz_sparse_mat_t A)
{
    slong pr, rk, pc, prev_pc, r;
    fmpz_t *ptr;

    /* Compute the rank (assuming HNF form) */
    for (rk = A->r; rk != 0; rk--)
    {
        if (!fmpz_sparse_vec_is_zero(&A->rows[rk - 1]) && LT(A, rk - 1).ind < A->c) break;
    }

    /* Check that first rk rows are in HNF */
    prev_pc = -1;
    for (pr = 0; pr < rk; pr++)
    {
        if (fmpz_sparse_vec_is_zero(&A->rows[pr])) return 0;
        pc = LT(A, pr).ind;
        if (pc >= A->c || pc <= prev_pc || fmpz_sgn(LT(A, pr).val) < 0) return 0;
        
        for (r = 0; r < pr; r++)
        {
            ptr = fmpz_sparse_vec_at(&A->rows[r], pc);
            if (ptr && ((fmpz_sgn(*ptr) < 0) || (fmpz_cmp(*ptr, LT(A, pr).val) >= 0))) return 0;
        }
        prev_pc = pc;
    }

    return 1;
}
