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

slong
nmod_mat_lu(slong * P, nmod_mat_t A, int rank_check)
{
    slong nrows, ncols, n, cutoff;
    int bits;
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

        const dot_params_t params = _nmod_vec_dot_params(n, A->mod);

        // TODO thresholds to re-examine after dot product changes
        if (params.method <= _DOT1                        // <= 0,1 limb
                || (params.method <= _DOT2 && n >= 12)    // <= 2 limbs (n >= 12 if exactly 2)
                || (params.method > _DOT2 && n >= 20))    // == 3 limbs && n >= 20
            return nmod_mat_lu_classical_delayed(P, A, rank_check);
        else
            return nmod_mat_lu_classical(P, A, rank_check);
    }
}
