/*
    Copyright (C) 2011, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"

/*
    Tuned for Zen 3 with BLAS. A gemm-backed nmod_mat_mul is now always
    available (nmod_mat_mul_blas uses flint_dgemm/flint_sgemm, which
    fall back to FLINT's own kernels when no external BLAS is
    configured), so the "without BLAS" table has been retired. It may be
    worth re-tuning against the fallback kernels.
*/
static const slong lu_cutoff_tab[64] = { 64, 64, 244, 280, 296, 320, 332,
    792, 800, 800, 800, 400, 400, 400, 400, 404, 404, 400, 408, 416, 400,
    412, 424, 408, 484, 1352, 1032, 260, 68, 56, 56, 56, 148, 160, 148, 156,
    156, 156, 120, 168, 160, 124, 156, 124, 156, 148, 112, 156, 148, 136,
    156, 160, 116, 136, 168, 148, 156, 160, 120, 128, 68, 64, 104, 104 };

slong
nmod_mat_lu(slong * P, nmod_mat_t A, int rank_check)
{
    slong nrows = A->r, ncols = A->c, n;
    n = FLINT_MIN(nrows, ncols);

    if (n <= 3 || (NMOD_BITS(A->mod) > 28 && n <= 7))
    {
        return nmod_mat_lu_classical(P, A, rank_check);
    }
    else
    {
        if (n < lu_cutoff_tab[NMOD_BITS(A->mod) - 1])
            return nmod_mat_lu_classical_delayed(P, A, rank_check);
        else
            return nmod_mat_lu_recursive(P, A, rank_check);
    }
}

