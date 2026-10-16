/*
    Copyright (C) 2023 Jean Kieffer
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

/*
    Compute the Smith normal form S of A and unimodular matrices U, V
    such that S = U * A * V.

    S is m x n, U is m x m (or NULL), V is n x n (or NULL).

    Algorithm: iterative Hermite normal form.
    Phase 1: alternately compute row-HNF and column-HNF (via transpose)
             until X is diagonal, tracking transforms in U and V.
    Phase 2: fix divisibility chain on the diagonal using xgcd,
             with O(m+n) row/column operations per pair instead of
             O((m+n)^3) matrix multiplications.

    Based on the implementation by Jean Kieffer in acb_theta, generalized
    to non-square matrices.
*/
void
fmpz_mat_snf_transform(fmpz_mat_t S, fmpz_mat_t U, fmpz_mat_t V,
    const fmpz_mat_t A)
{
    slong m = fmpz_mat_nrows(A);
    slong n = fmpz_mat_ncols(A);
    slong d = FLINT_MIN(m, n);
    slong max_iters, iter;
    fmpz_mat_t X, Xt, Mr, Mc;
    slong j, k, i;

    if (d == 0)
    {
        fmpz_mat_zero(S);
        if (U != NULL)
            fmpz_mat_one(U);
        if (V != NULL)
            fmpz_mat_one(V);
        return;
    }

    if (U == NULL && V == NULL)
    {
        fmpz_mat_snf(S, A);
        return;
    }

    fmpz_mat_init(X, m, n);
    fmpz_mat_init(Xt, n, m);
    fmpz_mat_init(Mr, m, m);
    fmpz_mat_init(Mc, n, n);

    fmpz_mat_set(X, A);
    if (U != NULL)
        fmpz_mat_one(U);
    if (V != NULL)
        fmpz_mat_one(V);

    max_iters = _fmpz_mat_snf_iter_bound(A);

    for (iter = 0; !fmpz_mat_is_diagonal(X); iter++)
    {
        if (iter >= max_iters)
            flint_throw(FLINT_ERROR,
                "(fmpz_mat_snf_transform): Phase 1 exceeded iteration "
                "bound (%wd).  Likely bug in fmpz_mat_hnf_transform or "
                "unexpected input; please report.\n", max_iters);

        /* Row HNF: X -> Mr*X, U -> Mr*U (Mr is m x m) */
        fmpz_mat_hnf_transform(X, Mr, X);
        if (U != NULL)
            fmpz_mat_mul(U, Mr, U);

        if (fmpz_mat_is_diagonal(X))
            break;

        /* Column HNF via transpose:
           Xt = X^T, compute Xt -> Mc*Xt (row HNF of transpose),
           then X = Xt^T, V -> V * Mc^T */
        fmpz_mat_transpose(Xt, X);
        fmpz_mat_hnf_transform(Xt, Mc, Xt);
        fmpz_mat_transpose(X, Xt);
        fmpz_mat_transpose(Mc, Mc);
        if (V != NULL)
            fmpz_mat_mul(V, V, Mc);
    }

    /*
        Phase 2: fix divisibility chain on diagonal.
        For each pair (j, k) with j < k, if X[j,j] does not divide X[k,k],
        apply xgcd-based row/column operations.

        Let a = X[j,j], b = X[k,k], d = gcd(a,b), a = d*p, b = d*q,
        u*a + v*b = d (so u*p + v*q = 1).

        Left transform L (rows j,k):
            L[j,j] = 1,  L[j,k] = v,  L[k,j] = q,  L[k,k] = vq-1
            det(L) = vq-1 - vq = -1

        Right transform R (cols j,k):
            R[j,j] = u,  R[j,k] = 1-up,  R[k,j] = 1,  R[k,k] = -p
            det(R) = -up - (1-up) = -1

        Result: L * diag(a,b) * R = diag(d, pqd) = diag(gcd, lcm).
    */
    {
        fmpz_t dd, uu, vv, pp, qq, vq_m1, one_m_up;
        fmpz * save_j, * save_k;

        fmpz_init(dd);
        fmpz_init(uu);
        fmpz_init(vv);
        fmpz_init(pp);
        fmpz_init(qq);
        fmpz_init(vq_m1);
        fmpz_init(one_m_up);

        save_j = _fmpz_vec_init(FLINT_MAX(m, n));
        save_k = _fmpz_vec_init(FLINT_MAX(m, n));

        for (j = 0; j < d; j++)
        {
            if (fmpz_is_one(fmpz_mat_entry(X, j, j)))
                continue;

            for (k = j + 1; k < d; k++)
            {
                if (fmpz_is_zero(fmpz_mat_entry(X, k, k)))
                    continue;
                if (!fmpz_is_zero(fmpz_mat_entry(X, j, j))
                    && fmpz_divisible(
                        fmpz_mat_entry(X, k, k),
                        fmpz_mat_entry(X, j, j)))
                    continue;

                fmpz_xgcd_canonical_bezout(dd, uu, vv,
                    fmpz_mat_entry(X, j, j),
                    fmpz_mat_entry(X, k, k));
                fmpz_divexact(pp,
                    fmpz_mat_entry(X, j, j), dd);
                fmpz_divexact(qq,
                    fmpz_mat_entry(X, k, k), dd);

                /* vq_m1 = v*q - 1 */
                fmpz_mul(vq_m1, vv, qq);
                fmpz_sub_ui(vq_m1, vq_m1, 1);

                /* one_m_up = 1 - u*p */
                fmpz_mul(one_m_up, uu, pp);
                fmpz_neg(one_m_up, one_m_up);
                fmpz_add_ui(one_m_up, one_m_up, 1);

                /*
                    Left transform on rows j,k of U:
                    new row j = old row j + v * old row k
                    new row k = q * old row j
                                + (vq-1) * old row k
                */
                if (U != NULL)
                {
                    fmpz * row_j = fmpz_mat_row(U, j);
                    fmpz * row_k = fmpz_mat_row(U, k);
                    _fmpz_vec_set(save_j, row_j, m);
                    _fmpz_vec_set(save_k, row_k, m);
                    _fmpz_vec_scalar_addmul_fmpz(row_j, save_k, m, vv);
                    _fmpz_vec_scalar_mul_fmpz(row_k, save_j, m, qq);
                    _fmpz_vec_scalar_addmul_fmpz(row_k, save_k, m, vq_m1);
                }

                /*
                    Right transform on cols j,k of V:
                    new col j = u * old col j + old col k
                    new col k = (1-up) * old col j
                                + (-p) * old col k
                */
                if (V != NULL)
                {
                    for (i = 0; i < n; i++)
                        fmpz_set(&save_j[i],
                            fmpz_mat_entry(V, i, j));
                    for (i = 0; i < n; i++)
                        fmpz_set(&save_k[i],
                            fmpz_mat_entry(V, i, k));

                    for (i = 0; i < n; i++)
                    {
                        fmpz_mul(
                            fmpz_mat_entry(V, i, j),
                            uu, &save_j[i]);
                        fmpz_add(
                            fmpz_mat_entry(V, i, j),
                            fmpz_mat_entry(V, i, j),
                            &save_k[i]);
                        fmpz_mul(
                            fmpz_mat_entry(V, i, k),
                            one_m_up, &save_j[i]);
                        fmpz_submul(
                            fmpz_mat_entry(V, i, k),
                            pp, &save_k[i]);
                    }
                }

                /* X[j,j] = gcd, X[k,k] = lcm */
                fmpz_set(fmpz_mat_entry(X, j, j), dd);
                fmpz_mul(fmpz_mat_entry(X, k, k),
                    pp, qq);
                fmpz_mul(fmpz_mat_entry(X, k, k),
                    fmpz_mat_entry(X, k, k), dd);
            }
        }

        _fmpz_vec_clear(save_k, FLINT_MAX(m, n));
        _fmpz_vec_clear(save_j, FLINT_MAX(m, n));
        fmpz_clear(dd);
        fmpz_clear(uu);
        fmpz_clear(vv);
        fmpz_clear(pp);
        fmpz_clear(qq);
        fmpz_clear(vq_m1);
        fmpz_clear(one_m_up);
    }

    /*
        Phase 3: fix signs.  SNF requires non-negative diagonal entries.
        If X[j,j] < 0, negate row j of U to absorb the sign.
    */
    for (j = 0; j < d; j++)
    {
        if (fmpz_sgn(fmpz_mat_entry(X, j, j)) < 0)
        {
            fmpz_neg(fmpz_mat_entry(X, j, j),
                fmpz_mat_entry(X, j, j));
            if (U != NULL)
            {
                fmpz * row_j = fmpz_mat_row(U, j);
                _fmpz_vec_neg(row_j, row_j, m);
            }
        }
    }

    fmpz_mat_set(S, X);

    fmpz_mat_clear(Mc);
    fmpz_mat_clear(Xt);
    fmpz_mat_clear(X);
    fmpz_mat_clear(Mr);
}
