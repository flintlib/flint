/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_MAT_IMPL_H
#define NMOD_MAT_IMPL_H

#include <stdio.h>
#include "thread_support.h"
#include "perm.h"
#include "nmod_mat.h"
#include "nmod_poly.h"

/*
with op = 0, computes D = A*B
with op = 1, computes D = C + A*B
with op = -1, computes D = C - A*B
*/

static __inline__ void
_nmod_mat_addmul_basic_op(mp_ptr * D, mp_ptr * const C, mp_ptr * const A,
    mp_ptr * const B, slong m, slong k, slong n, int op, nmod_t mod, int nlimbs)
{
    slong i, j;
    mp_limb_t c;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            c = _nmod_vec_dot_ptr(A[i], B, j, k, mod, nlimbs);

            if (op == 1)
                c = nmod_add(C[i][j], c, mod);
            else if (op == -1)
                c = nmod_sub(C[i][j], c, mod);

            D[i][j] = c;
        }
    }
}

static __inline__ int
_nmod_mat_pivot(nmod_mat_t A, slong start_row, slong col)
{
    slong j;
    mp_ptr u;

    if (nmod_mat_entry(A, start_row, col) != 0)
        return 1;

    for (j = start_row + 1; j < A->r; j++)
    {
        if (nmod_mat_entry(A, j, col) != 0)
        {
            u = A->rows[j];
            A->rows[j] = A->rows[start_row];
            A->rows[start_row] = u;

            return -1;
        }
    }
    return 0;
}

/* defined in det_howell.c */
int
_n_is_divisible(mp_ptr q, mp_limb_t b, mp_limb_t a, nmod_t N);

#endif
