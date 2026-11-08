/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_mat.h"

/*
    Applies the 2x2 matrix M (given in row-major order as a length-4
    vector) to rows i1, i2 of A, restricted to columns [c0, c1):

        [A[i1,k]]     [M[0]  M[1]] [A[i1,k]]
        [A[i2,k]]  <- [M[2]  M[3]] [A[i2,k]]

    i.e. A[[i1,i2], c0:c1] <- M * A[[i1,i2], c0:c1].

    Rows of a gr_mat are contiguous, so this is a unit-stride operation.
*/
int
gr_mat_apply_2x2_rows(gr_mat_t A, slong i1, slong i2, slong c0, slong c1, gr_srcptr M, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;

    if (c1 <= c0)
        return GR_SUCCESS;

    return _gr_vec_apply_2x2(GR_MAT_ENTRY(A, i1, c0, sz),
                             GR_MAT_ENTRY(A, i2, c0, sz),
                             1, c1 - c0, M, ctx);
}

/*
    Applies the transpose action of the 2x2 matrix M to columns j1, j2 of
    A, restricted to rows [r0, r1):

        A[r0:r1, [j1,j2]] <- A[r0:r1, [j1,j2]] * M^T

    that is, [A[k,j1], A[k,j2]] <- [M[0]*A[k,j1] + M[1]*A[k,j2],
                                    M[2]*A[k,j1] + M[3]*A[k,j2]].

    Columns are strided by A->stride.
*/
int
gr_mat_apply_2x2_cols(gr_mat_t A, slong j1, slong j2, slong r0, slong r1, gr_srcptr M, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;

    if (r1 <= r0)
        return GR_SUCCESS;

    return _gr_vec_apply_2x2(GR_MAT_ENTRY(A, r0, j1, sz),
                             GR_MAT_ENTRY(A, r0, j2, sz),
                             A->stride, r1 - r0, M, ctx);
}
