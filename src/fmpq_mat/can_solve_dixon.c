/*
    Copyright (C) 2022 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_mat.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq.h"
#include "fmpq_mat.h"
#include "fmpq_mat-impl.h"

/* Algorithm developed with Claus Fieker */

int
fmpq_mat_can_solve_fmpz_mat_dixon(fmpq_mat_t X,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    ulong p;
    fmpz_t tested;
    nmod_mat_t Ap, LU;
    int result = 0, success = 0;
    slong * perm, * pivots;
    slong i, j, k, col, rank;
    fmpz_mat_t Arank, Brank;
    fmpq_mat_t Xrank;
    fmpz_t det_bound;

    p = UWORD(1) << NMOD_MAT_OPTIMAL_MODULUS_BITS;
    fmpz_init(det_bound);
    fmpz_init(tested);
    fmpz_one(tested);
    nmod_mat_init(Ap, A->r, A->c, p);
    nmod_mat_init(LU, A->r, A->c, p);
    perm = flint_malloc(sizeof(slong)*A->r);
    pivots = flint_malloc(sizeof(slong)*A->c); /* only first rank entries are meaningful */

    fmpz_mat_det_bound(det_bound, A);

    while (1)
    {
        p = n_nextprime(p, 0);
        nmod_mat_set_mod(Ap, p);
        nmod_mat_set_mod(LU, p);

        fmpz_mat_get_nmod_mat(Ap, A);

        nmod_mat_set(LU, Ap);

        for (i = 0; i < A->r; i++)
            perm[i] = i;

        rank = nmod_mat_lu(perm, LU, 0);

        col = 0;
        for (i = 0; i < rank; i++)
        {
           while (nmod_mat_entry(LU, i, col) == 0)
              col++;

           pivots[i] = col;

           col++;
        }

        fmpz_mat_init(Arank, rank, rank);
        fmpz_mat_init(Brank, rank, B->c);
        fmpq_mat_init(Xrank, rank, B->c);

        for (i = 0; i < rank; i++)
        {
            k = 0;
            for (j = 0; j < A->c; j++)
            {
                if (k < rank && j == pivots[k])
                {
                    fmpz_set(fmpz_mat_entry(Arank, i, k), fmpz_mat_entry(A, perm[i], j));
                    k++;
                }
            }

            for (j = 0; j < B->c; j++)
                fmpz_set(fmpz_mat_entry(Brank, i, j), fmpz_mat_entry(B, perm[i], j));
        }

        success = fmpq_mat_solve_fmpz_mat_dixon(Xrank, Arank, Brank);

        if (success)
        {
            fmpq_mat_zero(X);

            for (i = 0; i < rank; i++)
            {
                for (j = 0; j < B->c; j++)
                    fmpq_set(fmpq_mat_entry(X, pivots[i], j), fmpq_mat_entry(Xrank, i, j));
            }

            result = _fmpq_mat_check_solution_fmpz_mat(X, A, B);
        }

        fmpz_mat_clear(Arank);
        fmpz_mat_clear(Brank);
        fmpq_mat_clear(Xrank);

        if (result)
            break;

        fmpz_mul_ui(tested, tested, p);
        if (fmpz_cmp(tested, det_bound) > 0)
            break;
    }

    fmpz_clear(det_bound);
    nmod_mat_clear(Ap);
    nmod_mat_clear(LU);
    fmpz_clear(tested);
    flint_free(perm);
    flint_free(pivots);

    return result;
}

int
fmpq_mat_can_solve_dixon(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
{
    fmpz_mat_t Anum;
    fmpz_mat_t Bnum;
    int success;

    if (A->r == 0 || B->c == 0)
    {
        fmpq_mat_zero(X);

        return 1;
    }

    if (A->c == 0)
    {
        fmpq_mat_zero(X);

        return fmpq_mat_is_zero(B);
    }

    fmpz_mat_init(Anum, A->r, A->c);
    fmpz_mat_init(Bnum, B->r, B->c);

    fmpq_mat_get_fmpz_mat_rowwise_2(Anum, Bnum, NULL, A, B);
    success = fmpq_mat_can_solve_fmpz_mat_dixon(X, Anum, Bnum);

    fmpz_mat_clear(Anum);
    fmpz_mat_clear(Bnum);

    return success;
}
