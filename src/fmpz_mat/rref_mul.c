/*
    Copyright (C) 2010-2012, 2026 Fredrik Johansson
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "perm.h"
#include "ulong_extras.h"
#include "nmod_mat.h"
#include "fmpz.h"
#include "fmpz_mat.h"

int
fmpz_mat_rref_upper_certify_lu_mod_p(fmpz_mat_t E, fmpz_t den, const fmpz_mat_t A,
                          slong rank, const slong * P, const slong * pivs)
{
    slong i, j, m, n;
    fmpz_mat_t B, C, D, E2, F, FD;

    m = fmpz_mat_nrows(A);
    n = fmpz_mat_ncols(A);

    fmpz_mat_init(B, rank, rank);
    fmpz_mat_init(C, rank, n - rank);

    /* B = pivot rows x pivot columns, C = pivot rows x non-pivot columns */
    for (i = 0; i < rank; i++)
    {
        for (j = 0; j < rank; j++)
            fmpz_set(fmpz_mat_entry(B, i, j),
                    fmpz_mat_entry(A, P[i], pivs[j]));
        for (j = 0; j < n - rank; j++)
            fmpz_set(fmpz_mat_entry(C, i, j),
                    fmpz_mat_entry(A, P[i], pivs[rank + j]));
    }

    /* solve B*E2 = den*C (B is invertible mod p, hence over Z) */
    fmpz_mat_init(E2, rank, n - rank);
    if (!fmpz_mat_solve(E2, den, B, C))
    {
        flint_throw(FLINT_ERROR, "(fmpz_mat_rref_upper_certify_lu_mod_p): "
                "Singular input matrix for solve.\n");
    }

    fmpz_mat_clear(B);
    fmpz_mat_clear(C);

    /* assemble the scaled rref rows into E: den on the pivot "diagonal", E2 in
       the non-pivot columns; the remaining pivot-column entries stay zero */
    for (i = 0; i < rank; i++)
    {
        fmpz_set(fmpz_mat_entry(E, i, pivs[i]), den);
        for (j = 0; j < n - rank; j++)
            fmpz_set(fmpz_mat_entry(E, i, pivs[rank + j]),
                    fmpz_mat_entry(E2, i, j));
    }
    fmpz_mat_clear(E2);

    if (!fmpz_mat_is_in_rref_with_rank(E, den, rank))
        return 0;

    /* D = nullspace basis for E (n x (n - rank)) */
    fmpz_mat_init(D, n, n - rank);
    for (j = 0; j < n - rank; j++)
    {
        fmpz_set(fmpz_mat_entry(D, pivs[rank + j], j), den);
        for (i = 0; i < rank; i++)
            fmpz_neg(fmpz_mat_entry(D, pivs[i], j),
                    fmpz_mat_entry(E, i, pivs[rank + j]));
    }

    /* F = the non-pivot rows of A. The rank is certified iff F*D == 0, i.e.
       every remaining row lies in the span of the pivot rows. */
    fmpz_mat_init(F, m - rank, n);
    for (i = 0; i < m - rank; i++)
        for (j = 0; j < n; j++)
            fmpz_set(fmpz_mat_entry(F, i, j),
                    fmpz_mat_entry(A, P[rank + i], j));

    fmpz_mat_init(FD, m - rank, n - rank);
    fmpz_mat_mul(FD, F, D);
    fmpz_mat_clear(F);
    fmpz_mat_clear(D);

    if (!fmpz_mat_is_zero(FD))
    {
        fmpz_mat_clear(FD);
        return 0;
    }

    fmpz_mat_clear(FD);
    return 1;
}

int
fmpz_mat_rank_certify_lu_mod_p(const fmpz_mat_t A,
    slong rank, const slong * P, const slong * pivs)
{
    int result;

    fmpz_mat_t E;
    fmpz_t den;

    fmpz_init(den);
    fmpz_mat_init(E, rank, A->c);

    result = fmpz_mat_rref_upper_certify_lu_mod_p(E, den, A, rank, P, pivs);

    fmpz_clear(den);
    fmpz_mat_clear(E);

    return result;
}

slong
fmpz_mat_rref_mul(fmpz_mat_t R, fmpz_t den, const fmpz_mat_t A)
{
    nmod_mat_t Amod;
    fmpz_mat_t E;
    ulong p;
    slong i, j, m, n, rank, * pivs, * P;

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

        rank = nmod_mat_lu_with_pivots(P, pivs, Amod);
        nmod_mat_clear(Amod);

        /* stop early if the rank is the number of columns: the rref is the
           identity (with zero rows below) and the denominator is 1 */
        if (rank == n)
        {
            fmpz_mat_one(R);
            fmpz_one(den);
            break;
        }

        /* otherwise compute and certify the rref for this prime */
        fmpz_mat_init(E, rank, n);
        if (fmpz_mat_rref_upper_certify_lu_mod_p(E, den, A, rank, P, pivs))
        {
            /* E is valid */
            /* write the entries of E into R and zeroes at the bottom */
            for (i = 0; i < rank; i++)
                for (j = 0; j < n; j++)
                    fmpz_set(fmpz_mat_entry(R, i, j), fmpz_mat_entry(E, i, j));
            for (i = rank; i < m; i++)
                for (j = 0; j < n; j++)
                    fmpz_zero(fmpz_mat_entry(R, i, j));
            fmpz_mat_clear(E);
            break;
        }
        else
        {
            fmpz_mat_clear(E);          /* unlucky prime: discard and try again */
        }
    }

    flint_free(pivs);
    _perm_clear(P);

    return rank;
}
