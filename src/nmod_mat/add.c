/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_mat.h"

void
nmod_mat_add(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    slong i;

    FLINT_ASSERT(C->r == A->r);
    FLINT_ASSERT(C->r == B->r);
    FLINT_ASSERT(C->c == A->c);
    FLINT_ASSERT(C->c == B->c);

    if (C->c == 0)
        return;

    for (i = 0; i < C->r; i++)
    {
        _nmod_vec_add(nmod_mat_entry_ptr(C, i, 0),
            nmod_mat_entry_ptr(A, i, 0), nmod_mat_entry_ptr(B, i, 0), C->c, C->mod);
    }
}
