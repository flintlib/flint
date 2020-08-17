/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2020 William Hart

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
#include "thread_support.h"

void
nmod_mat_mul(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    slong m, k, n, cutoff;

    m = A->r;
    k = A->c;
    n = B->c;

    if (C == A || C == B)
    {
        nmod_mat_t T;
        nmod_mat_init(T, m, n, A->mod.n);
        nmod_mat_mul(T, A, B);
        nmod_mat_swap(C, T);
        nmod_mat_clear(T);
        return;
    }

    if (FLINT_BITS == 64 && C->mod.n < 2048)
        cutoff = 400;
    else
        cutoff = 200;

    if (flint_get_num_threads() > 1)
	nmod_mat_mul_classical_threaded(C, A, B);
    else if (m < cutoff || n < cutoff || k < cutoff)
        nmod_mat_mul_classical(C, A, B);
    else
        nmod_mat_mul_strassen(C, A, B);
}
