/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

slong
TEMPLATE(T, mat_nullspace) (TEMPLATE(T, mat_t) X, const TEMPLATE(T, mat_t) A,
                            const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j, k, m, n, rank, nullity;
    slong *p;
    slong *pivots;
    slong *nonpivots;
    TEMPLATE(T, mat_t) tmp;

    m = A->r;
    n = A->c;

    p = flint_malloc(sizeof(slong) * FLINT_MAX(m, n));

    TEMPLATE(T, mat_init_set) (tmp, A, ctx);
    rank = TEMPLATE(T, mat_rref) (tmp, ctx);
    nullity = n - rank;

    TEMPLATE(T, mat_zero) (X, ctx);

    if (rank == 0)
    {
        for (i = 0; i < nullity; i++)
            TEMPLATE(T, one) (TEMPLATE(T, mat_entry) (X, i, i), ctx);
    }
    else if (nullity)
    {
        pivots = p;             /* length = rank */
        nonpivots = p + rank;   /* length = nullity */

        for (i = j = k = 0; i < rank; i++)
        {
            while (TEMPLATE(T, is_zero)
                   (TEMPLATE(T, mat_entry) (tmp, i, j), ctx))
            {
                nonpivots[k] = j;
                k++;
                j++;
            }
            pivots[i] = j;
            j++;
        }
        while (k < nullity)
        {
            nonpivots[k] = j;
            k++;
            j++;
        }

        for (i = 0; i < nullity; i++)
        {
            for (j = 0; j < rank; j++)
            {
                TEMPLATE(T, neg) (TEMPLATE(T, mat_entry) (X, pivots[j], i),
                                  TEMPLATE(T, mat_entry) (tmp, j,
                                                          nonpivots[i]), ctx);
            }

            TEMPLATE(T, one) (TEMPLATE(T, mat_entry) (X, nonpivots[i], i),
                              ctx);
        }
    }

    flint_free(p);
    TEMPLATE(T, mat_clear) (tmp, ctx);

    return nullity;
}

#endif
