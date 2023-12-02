/*
    Copyright (C) 2010-2012 Fredrik Johansson
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "perm.h"
#include "nmod_mat.h"
#include "fmpz.h"
#include "fmpz_mat.h"

slong
fmpz_mat_rref_mul(fmpz_mat_t R, fmpz_t den, const fmpz_mat_t A)
{
    nmod_mat_t Amod;
    mp_limb_t p;
    slong i, j, m, n, rank, * pivs, * P;
    fmpz_mat_t B, C, D, E, E2, F, FD;

    m = fmpz_mat_nrows(A);
    n = fmpz_mat_ncols(A);

    pivs = (slong *) flint_malloc(n * sizeof(slong));
    P = _perm_init(m);

    /* use 16 bit primes to ensure it is unlikely we hit a bad one and so that
       the modular computations are not too long */
    p = 1 << 16;

    while (1)
    {
        p = n_nextprime(p, 1);
        nmod_mat_init(Amod, m, n, p);
        fmpz_mat_get_nmod_mat(Amod, A);

        rank = _nmod_mat_rref(Amod, pivs, P);

        nmod_mat_clear(Amod);

        /* stop early if the rank is the number of columns */
        if (rank == n)
        {
            fmpz_mat_one(R);
            fmpz_one(den);

            flint_free(pivs);
            _perm_clear(P);

            return rank;
        }

        fmpz_mat_init(B, rank, rank);
        fmpz_mat_init(C, rank, n - rank);

        /* set B to be the pivot columns and rows and C to be the non-pivot
           columns in the pivot rows */
        for (i = 0; i < rank; i++)
        {
            for (j = 0; j < rank; j++)
                fmpz_set(fmpz_mat_entry(B, i, j),
                        fmpz_mat_entry(A, P[i], pivs[j]));
            for (j = 0; j < n - rank; j++)
                fmpz_set(fmpz_mat_entry(C, i, j),
                        fmpz_mat_entry(A, P[i], pivs[rank + j]));
        }

        /* solve B*E2 = den*C */
        fmpz_mat_init(E2, rank, n - rank);
        if (!fmpz_mat_solve(E2, den, B, C))
        {
            flint_throw(FLINT_ERROR, "(fmpz_mat_rref_mul): "
                    "Singular input matrix for solve.\n");
        }

        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_init(E, rank, n);

        /* move columns of E2 and identity matrix into E so that it should be
           in rref */
        for (i = 0; i < rank; i++)
        {
            fmpz_set(fmpz_mat_entry(E, i, pivs[i]), den);
            for (j = 0; j < n - rank; j++)
                fmpz_set(fmpz_mat_entry(E, i, pivs[rank + j]),
                        fmpz_mat_entry(E2, i, j));
        }
        fmpz_mat_clear(E2);

        if (!fmpz_mat_is_in_rref_with_rank(E, den, rank))
        {
            fmpz_mat_clear(E);
            continue;
        }

        /* set D to be the nullspace basis matrix for E */
        fmpz_mat_init(D, n, n - rank);

        for (j = 0; j < n - rank; j++)
        {
            fmpz_set(fmpz_mat_entry(D, pivs[rank + j], j), den);
            for (i = 0; i < rank; i++)
                fmpz_neg(fmpz_mat_entry(D, pivs[i], j),
                        fmpz_mat_entry(E, i, pivs[rank + j]));
        }

        fmpz_mat_init(F, m - rank, n);

        for (i = 0; i < m - rank; i++)
            for (j = 0; j < n; j++)
                fmpz_set(fmpz_mat_entry(F, i, j),
                        fmpz_mat_entry(A, P[rank + i], j));

        fmpz_mat_init(FD, m - rank, n - rank);
        fmpz_mat_mul(FD, F, D);
        fmpz_mat_clear(F);
        fmpz_mat_clear(D);

        /* if FD = 0 we have computed the rref right so stop, otherwise try a
           different p in the next iteration */
        if (fmpz_mat_is_zero(FD))
            break;

        fmpz_mat_clear(E);
        fmpz_mat_clear(FD);
    }

    /* write the entries of E into R and zeroes at the bottom */
    for (i = 0; i < rank; i++)
        for (j = 0; j < n; j++)
            fmpz_set(fmpz_mat_entry(R, i, j), fmpz_mat_entry(E, i, j));
    for (i = rank; i < m; i++)
        for (j = 0; j < n; j++)
            fmpz_zero(fmpz_mat_entry(R, i, j));

    fmpz_mat_clear(E);
    fmpz_mat_clear(FD);
    flint_free(pivs);
    _perm_clear(P);

    return rank;
}

