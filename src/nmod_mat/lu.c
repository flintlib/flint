/*
    Copyright (C) 2011, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_mat.h"

slong 
nmod_mat_lu(slong * P, nmod_mat_t A, int rank_check)
{
    slong nrows, ncols, n, cutoff;
    int nlimbs, bits;
    nrows = A->r;
    ncols = A->c;

    n = FLINT_MIN(nrows, ncols);

    if (n <= 3)
    {
        return nmod_mat_lu_classical(P, A, rank_check);
    }
    else
    {
        if (n >= 20)
        {
            bits = NMOD_BITS(A->mod);

            if (bits >= FLINT_BITS - 1)
                cutoff = 80;
            else if (bits >= FLINT_BITS / 2 - 2)
                cutoff = 60;
            else if (bits >= FLINT_BITS / 4 - 1)
                cutoff = 180;
            else
                cutoff = 60;

            if (n >= cutoff)
                return nmod_mat_lu_recursive(P, A, rank_check);
        }

        nlimbs = _nmod_vec_dot_bound_limbs(n, A->mod);

        if (nlimbs <= 1 || (nlimbs == 2 && n >= 12) || (nlimbs == 3 && n >= 20))
            return nmod_mat_lu_classical_delayed(P, A, rank_check);
        else
            return nmod_mat_lu_classical(P, A, rank_check);
    }
}
