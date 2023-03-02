/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"
#include "nmod_vec.h"

void
nmod_mat_submul(nmod_mat_t D, const nmod_mat_t C,
                                const nmod_mat_t A, const nmod_mat_t B)
{
    slong m, k, n, cutoff;

    m = A->r;
    k = A->c;
    n = B->c;

    if (FLINT_BITS == 64 && C->mod.n < 2048)
        cutoff = 400;
    else
        cutoff = 200;

    if (flint_get_num_threads() == 1 && (m < cutoff || n < cutoff || k < cutoff))
    {
        _nmod_mat_mul_classical_op(D, C, A, B, -1);
    }
    else
    {
        nmod_mat_t tmp;
        nmod_mat_init(tmp, m, n, A->mod.n);
        nmod_mat_mul(tmp, A, B);
        nmod_mat_sub(D, C, tmp);
        nmod_mat_clear(tmp);
    }
}
