/*
    Copyright (C) 2014 Alex J. Best
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"

/*
    Compute SNF via iterative Hermite normal form.

    Phase 1: alternately compute row-HNF and column-HNF (via transpose)
             until the matrix is diagonal.
    Phase 2: fix divisibility chain on the diagonal using gcd/lcm.
    Phase 3: make diagonal entries non-negative.

    This is the same algorithm as fmpz_mat_snf_transform but without
    tracking the transformation matrices U and V.
*/
static void
_fmpz_mat_snf_iterative_hermite(fmpz_mat_t S, const fmpz_mat_t A)
{
    slong m = fmpz_mat_nrows(A);
    slong n = fmpz_mat_ncols(A);
    slong d = FLINT_MIN(m, n);
    slong j, k;
    slong max_iters, iter;
    fmpz_mat_t X, Xt;
    fmpz_t dd, pp, qq;

    if (d == 0)
    {
        fmpz_mat_zero(S);
        return;
    }

    fmpz_mat_init(X, m, n);
    fmpz_mat_init(Xt, n, m);
    fmpz_init(dd);
    fmpz_init(pp);
    fmpz_init(qq);

    fmpz_mat_set(X, A);

    max_iters = _fmpz_mat_snf_iter_bound(A);

    for (iter = 0; !fmpz_mat_is_diagonal(X); iter++)
    {
        if (iter >= max_iters)
            flint_throw(FLINT_ERROR,
                "(fmpz_mat_snf): Phase 1 exceeded iteration bound "
                "(%wd).  Likely bug in fmpz_mat_hnf or unexpected "
                "input; please report.\n", max_iters);

        /* Row HNF */
        fmpz_mat_hnf(X, X);

        if (fmpz_mat_is_diagonal(X))
            break;

        /* Column HNF via transpose */
        fmpz_mat_transpose(Xt, X);
        fmpz_mat_hnf(Xt, Xt);
        fmpz_mat_transpose(X, Xt);
    }

    /*
        Phase 2: fix divisibility chain on diagonal.
        For each pair (j, k) with j < k, if X[j,j] does not divide X[k,k],
        replace with gcd and lcm.
    */
    for (j = 0; j < d; j++)
    {
        if (fmpz_is_one(fmpz_mat_entry(X, j, j)))
            continue;

        for (k = j + 1; k < d; k++)
        {
            if (fmpz_is_zero(fmpz_mat_entry(X, k, k)))
                continue;
            if (!fmpz_is_zero(fmpz_mat_entry(X, j, j))
                && fmpz_divisible(fmpz_mat_entry(X, k, k),
                                  fmpz_mat_entry(X, j, j)))
                continue;

            fmpz_gcd(dd, fmpz_mat_entry(X, j, j),
                fmpz_mat_entry(X, k, k));
            fmpz_divexact(pp, fmpz_mat_entry(X, j, j), dd);
            fmpz_divexact(qq, fmpz_mat_entry(X, k, k), dd);

            /* X[j,j] = gcd, X[k,k] = lcm = p*q*d */
            fmpz_set(fmpz_mat_entry(X, j, j), dd);
            fmpz_mul(fmpz_mat_entry(X, k, k), pp, qq);
            fmpz_mul(fmpz_mat_entry(X, k, k),
                fmpz_mat_entry(X, k, k), dd);
        }
    }

    /*
        Phase 3: make diagonal entries non-negative.
    */
    for (j = 0; j < d; j++)
    {
        if (fmpz_sgn(fmpz_mat_entry(X, j, j)) < 0)
            fmpz_neg(fmpz_mat_entry(X, j, j),
                fmpz_mat_entry(X, j, j));
    }

    fmpz_mat_set(S, X);

    fmpz_mat_clear(Xt);
    fmpz_mat_clear(X);
    fmpz_clear(dd);
    fmpz_clear(pp);
    fmpz_clear(qq);
}

void
fmpz_mat_snf(fmpz_mat_t S, const fmpz_mat_t A)
{
    if (fmpz_mat_is_diagonal(A))
    {
        fmpz_mat_snf_diagonal(S, A);
        return;
    }

    _fmpz_mat_snf_iterative_hermite(S, A);
}
