/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include "perm.h"

slong
TEMPLATE(T, mat_rref) (TEMPLATE(T, mat_t) B, const TEMPLATE(T, mat_t) A, const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j, k, n, rank;
    slong *pivots;
    slong *nonpivots;
    slong *P;
    TEMPLATE(T, struct) * e;
    TEMPLATE(T, mat_t) U, V;

    TEMPLATE(T, mat_set)(B, A, ctx);

    if (TEMPLATE(T, mat_is_zero)(B, ctx))
        return 0;

    if (B->r == 1)
    {
        TEMPLATE(T, struct) * c;
        slong r = 0;

        for (i = 0; i < B->c; i++)
        {
            c = TEMPLATE(T, mat_entry)(B, 0, i);
            if (!TEMPLATE(T, is_zero)(c, ctx))
            {
                r = 1;
                if (TEMPLATE(T, is_one)(c, ctx))
                    break;

                TEMPLATE(T, inv)(c, c, ctx);
                for (j = i + 1;j < B->c; j++)
                {
                    TEMPLATE(T, mul)(TEMPLATE(T, mat_entry)(B, 0, j), TEMPLATE(T, mat_entry)(B, 0, j), c, ctx);
                }
                TEMPLATE(T, one)(c, ctx);
                break;
            }
        }
        return r;
    }

    n = A->c;

    P = _perm_init(TEMPLATE(T, mat_nrows) (B, ctx));
    rank = TEMPLATE(T, mat_lu) (P, B, 0, ctx);
    _perm_clear(P);

    if (rank == 0)
        return rank;

    /* Clear L */
    for (i = 0; i < B->r; i++)
        for (j = 0; j < FLINT_MIN(i, rank); j++)
            TEMPLATE(T, zero) (TEMPLATE(T, mat_entry) (B, i, j), ctx);

    /* We now reorder U to proper upper triangular form U | V
       with U full-rank triangular, set V = U^(-1) V, and then
       put the column back in the original order.

       An improvement for some matrices would be to compress V by
       discarding columns containing nothing but zeros. */

    TEMPLATE(T, mat_init) (U, rank, rank, ctx);
    TEMPLATE(T, mat_init) (V, rank, n - rank, ctx);

    pivots = flint_malloc(sizeof(slong) * rank);
    nonpivots = flint_malloc(sizeof(slong) * (n - rank));

    for (i = j = k = 0; i < rank; i++)
    {
        while (TEMPLATE(T, is_zero) (TEMPLATE(T, mat_entry) (B, i, j), ctx))
        {
            nonpivots[k] = j;
            k++;
            j++;
        }
        pivots[i] = j;
        j++;
    }
    while (k < n - rank)
    {
        nonpivots[k] = j;
        k++;
        j++;
    }

    for (i = 0; i < rank; i++)
    {
        for (j = 0; j <= i; j++)
        {
            e = TEMPLATE(T, mat_entry) (B, j, pivots[i]);
            TEMPLATE(T, mat_entry_set) (U, j, i, e, ctx);
        }
    }

    for (i = 0; i < n - rank; i++)
    {
        for (j = 0; j < rank; j++)
        {
            e = TEMPLATE(T, mat_entry) (B, j, nonpivots[i]);
            TEMPLATE(T, mat_entry_set) (V, j, i, e, ctx);
        }
    }

    TEMPLATE(T, mat_solve_triu) (V, U, V, 0, ctx);

    /* Clear pivot columns */
    for (i = 0; i < rank; i++)
    {
        for (j = 0; j <= i; j++)
        {
            if (i == j)
            {
                TEMPLATE(T, one) (TEMPLATE(T, mat_entry) (B, j, pivots[i]),
                                  ctx);
            }
            else
            {
                TEMPLATE(T, zero) (TEMPLATE(T, mat_entry) (B, j, pivots[i]),
                                   ctx);
            }
        }
    }

    /* Write back the actual content */
    for (i = 0; i < n - rank; i++)
    {
        for (j = 0; j < rank; j++)
            TEMPLATE(T, mat_entry_set) (B, j, nonpivots[i],
                                        TEMPLATE(T, mat_entry) (V, j, i), ctx);
    }

    TEMPLATE(T, mat_clear) (U, ctx);
    TEMPLATE(T, mat_clear) (V, ctx);

    flint_free(pivots);
    flint_free(nonpivots);

    return rank;
}


#endif
