/*
    Copyright 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

/*
    Follows section 13.4 of Siegfried M. Rump, "Verication methods:
    Rigorous results using floating-point arithmetic",
    Acta Numerica 19 (2010), pp. 287 - 449, implemented as
    verifyeig() in INTLAB.

    Cheat sheet for the formulas in Rump's paper.
    Assuming U = first n-k indices, V = last k indices.

    U U^T M   = selects first n-k rows from M,    zeroing rest
    V V^T M   = selects last k rows from M,       zeroing rest
    M U U^T   = selects first n-k columns from M, zeroing rest
    M V V^T   = selects last k columns from M,    zeroing rest

    U^T M     = selects first n-k rows from M, truncating matrix
    V^T M     = selects last k rows from M, truncating matrix
    M U       = selects first n-k columns from M, truncating matrix
    M V       = selects last k columns from M, truncating matrix

    X V^T     = extends X to n x n matrix (placing X on right)
*/

static void
acb_approx_neg(acb_t res, const acb_t x)
{
    arf_neg(arb_midref(acb_realref(res)), arb_midref(acb_realref(x)));
    arf_neg(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(x)));
}

static void
acb_approx_sub(acb_t res, const acb_t x, const acb_t y, slong prec)
{
    arf_sub(arb_midref(acb_realref(res)), arb_midref(acb_realref(x)), arb_midref(acb_realref(y)), prec, ARF_RND_DOWN);
    arf_sub(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(x)), arb_midref(acb_imagref(y)), prec, ARF_RND_DOWN);
}

/* todo: separate out */
void
acb_mat_bound_max_norm(mag_t res, const acb_mat_t A)
{
    mag_t t;
    slong i, j;

    mag_init(t);
    mag_zero(res);

    for (i = 0; i < acb_mat_nrows(A); i++)
    {
        for (j = 0; j < acb_mat_ncols(A); j++)
        {
            acb_get_mag(t, acb_mat_entry(A, i, j));
            mag_max(res, res, t);
        }
    }

    mag_clear(t);
}

static void
arb_mat_nonnegative_eig_bound(mag_t eps, const arb_mat_t M, slong prec)
{
    /* Cheap, but poor for defective eigenvalues */
    arb_mat_bound_frobenius_norm(eps, M);

    /* Use Perron root bound. TODO: do something direct for k = 2. */
    if (1)
    {
        acb_mat_t A, R, E;
        arb_mat_t V, MV;
        mag_t tm, um, vbound;
        slong i, j, k;

        k = arb_mat_nrows(M);

        acb_mat_init(A, k, k);
        acb_mat_init(R, k, k);
        acb_mat_init(E, 1, k);
        arb_mat_init(V, k, k);
        arb_mat_init(MV, k, k);
        mag_init(tm);
        mag_init(um);
        mag_init(vbound);

        acb_mat_set_arb_mat(A, M);
        /* TODO: could probably lower precision if precision is very high? */
        acb_mat_approx_eig_qr(acb_mat_entry(E, 0, 0), NULL, R, A, NULL, 0, prec);

        for (i = 0; i < k; i++)
        {
            for (j = 0; j < k; j++)
            {
                acb_get_mag(tm, acb_mat_entry(R, i, j));
                arf_set_mag(arb_midref(arb_mat_entry(V, i, j)), tm);
            }
        }

        arb_mat_mul(MV, M, V, MAG_BITS);

        for (j = 0; j < k; j++)
        {
            mag_zero(vbound);

            for (i = 0; i < k; i++)
            {
                arb_get_mag(tm, arb_mat_entry(MV, i, j));
                arb_get_mag_lower(um, arb_mat_entry(V, i, j));
                mag_div(tm, tm, um);
                mag_max(vbound, vbound, tm);
            }

            mag_min(eps, eps, vbound);
        }

        acb_mat_clear(A);
        acb_mat_clear(R);
        acb_mat_clear(E);
        arb_mat_clear(V);
        arb_mat_clear(MV);
        mag_clear(tm);
        mag_clear(um);
        mag_clear(vbound);
    }
}

static void
acb_approx_mag(mag_t res, const acb_t x)
{
    mag_t t;
    mag_init(t);
    arf_get_mag(res, arb_midref(acb_realref(x)));
    arf_get_mag(t, arb_midref(acb_imagref(x)));
    mag_hypot(res, res, t);
    mag_clear(t);
}

/* Extract k largest rows to freeze */
static void
partition_X_sorted(slong * u, slong * v, const acb_mat_t X, slong prec)
{
    slong i, j, n, k, c;
    slong * row_idx;
    mag_ptr row_mag;
    mag_t t;

    n = acb_mat_nrows(X);
    k = acb_mat_ncols(X);

    row_mag = _mag_vec_init(n);
    row_idx = flint_malloc(sizeof(slong) * n);
    mag_init(t);

    for (i = 0; i < n; i++)
    {
        row_idx[i] = i;

        for (j = 0; j < k; j++)
        {
            acb_approx_mag(t, acb_mat_entry(X, i, j));
            mag_add(row_mag + i, row_mag + i, t);
        }
    }

    /* Bubble sort... */
    for (i = 0; i < n - 1; i++)
    {
        for (j = 0; j < n - i - 1; j++)
        {
            if (mag_cmp(row_mag + j, row_mag + j + 1) > 0)
            {
                mag_swap(row_mag + j, row_mag + j + 1);
                c = row_idx[j];
                row_idx[j] = row_idx[j + 1];
                row_idx[j + 1] = c;
            }
        }
    }

    /* Not frozen rows of the approximation. */
    for (i = 0; i < n - k; i++)
        u[i] = row_idx[i];

    /* Frozen rows of the approximation. */
    for (i = 0; i < k; i++)
        v[i] = row_idx[n - k + i];

    _mag_vec_clear(row_mag, n);
    flint_free(row_idx);
    mag_clear(t);
}

static void
partition_X_trivial(slong * u, slong * v, const acb_mat_t X, slong prec)
{
    slong n, k, i;

    n = acb_mat_nrows(X);
    k = acb_mat_ncols(X);

    /* Not frozen rows of the approximation. */
    for (i = 0; i < n - k; i++)
        u[i] = i;
    /* Frozen rows of the approximation. */
    for (i = 0; i < k; i++)
        v[i] = n - k + i;
}

void
acb_mat_eig_enclosure_rump(acb_t lambda, acb_mat_t J, acb_mat_t X, const acb_mat_t A,
    const acb_t lambda_approx, const acb_mat_t X_approx, slong prec)
{
    slong n, k, i, j, iter, maxiter;
    slong *u, *v;
    acb_mat_t R, I, T, Y, Y0, UY, VY, Yeps;
    mag_t eps;

    n = acb_mat_nrows(A);
    k = acb_mat_ncols(X_approx);

    if (k < 1 || k > n || n != acb_mat_nrows(X_approx) || n != acb_mat_ncols(A))
    {
        flint_throw(FLINT_ERROR, "bad matrix dimensions in acb_mat_eig_enclosure_rump\n");
    }

    /* Not frozen rows of the approximation. */
    u = flint_malloc(sizeof(slong) * (n - k));
    /* Frozen rows of the approximation. */
    v = flint_malloc(sizeof(slong) * k);

    if (k == 1)
        partition_X_sorted(u, v, X_approx, prec);
    else
        partition_X_trivial(u, v, X_approx, prec);

    mag_init(eps);
    acb_mat_init(R, n, n);
    acb_mat_init(UY, n, k);
    acb_mat_init(VY, k, k);
    acb_mat_init(T, n, n);
    acb_mat_init(Y, n, k);
    acb_mat_init(Y0, n, k);
    acb_mat_init(Yeps, n, k);

    /* Preconditioner:
       R ~= ((A - lambda_approx I) U U^T - X_approx V^T)^(-1) */
    acb_mat_get_mid(R, A);
    for (i = 0; i < n; i++)
        acb_approx_sub(acb_mat_entry(R, i, i),
            acb_mat_entry(R, i, i), lambda_approx, prec);

    for (i = 0; i < n; i++)
        for (j = 0; j < k; j++)
            acb_approx_neg(acb_mat_entry(R, i, v[j]),
                acb_mat_entry(X_approx, i, j));

    acb_mat_init(I, n, n);
    acb_mat_one(I);
    acb_mat_approx_solve(R, R, I, prec);
    acb_mat_clear(I);

    /* T = I - R * ((A - lambda_approx I) U U^T - X_approx V^T) */
    /* Y = Y_0 = -R * ((A - lambda_approx I) X_approx) */
    acb_mat_set(T, A);
    for (i = 0; i < n; i++)
        acb_sub(acb_mat_entry(T, i, i), acb_mat_entry(T, i, i), lambda_approx, prec);

    acb_mat_mul(Y0, T, X_approx, prec);
    acb_mat_mul(Y0, R, Y0, prec);
    acb_mat_neg(Y0, Y0);
    acb_mat_set(Y, Y0);

    for (i = 0; i < n; i++)
        for (j = 0; j < k; j++)
            acb_neg(acb_mat_entry(T, i, v[j]), acb_mat_entry(X_approx, i, j));

    acb_mat_mul(T, R, T, prec);
    acb_mat_neg(T, T);
    for (i = 0; i < n; i++)
        acb_add_ui(acb_mat_entry(T, i, i), acb_mat_entry(T, i, i), 1, prec);

    /* Iteration with epsilon-inflation */
    /* Y represents the error with respect to lambda_approx and X_approx */
    /* TODO: what number of iterations is actually reasonable? */
    /* TODO: what size of epsilon is actually reasonable? */
    maxiter = 5 + FLINT_BIT_COUNT(prec);
    for (iter = 0; iter < maxiter; iter++)
    {
        /* Inflate Y. TODO: make it elementwise? */
        acb_mat_bound_max_norm(eps, Y);
        if (mag_is_zero(eps))
            mag_set_ui_2exp_si(eps, 1, -20 * prec);
        mag_mul_2exp_si(eps, eps, -3 + 2 * iter);
        /* if (iter > 3)
            mag_mul_2exp_si(eps, eps, (prec / 2) * (iter - 3) / (maxiter - 3)); */

        acb_mat_add_error_mag(Y, eps);
        acb_mat_set(Yeps, Y);

        /* Y = Y0 + T Y + R ((U U^T Y) V^T Y) */
        acb_mat_zero(UY);
        acb_mat_zero(VY);

        /* U U^T Y -- zero the rows at indices v. */
        acb_mat_set(UY, Y);
        for (i = 0; i < k; i++)
            for (j = 0; j < k; j++)
                acb_zero(acb_mat_entry(UY, v[i], j));

        /* V^T Y -- extract rows at indices v */
        for (i = 0; i < k; i++)
            for (j = 0; j < k; j++)
                acb_set(acb_mat_entry(VY, i, j), acb_mat_entry(Y, v[i], j));

        acb_mat_mul(UY, UY, VY, prec);
        acb_mat_mul(UY, R, UY, prec);

        acb_mat_mul(Y, T, Y, prec);
        acb_mat_add(Y, Y, UY, prec);
        acb_mat_add(Y, Y, Y0, prec);

        if (acb_mat_contains(Yeps, Y))
        {
            acb_get_mid(lambda, lambda_approx);

            if (J != NULL)
            {
                /* J = lambda_approx I_k + V^T Y */
                for (i = 0; i < k; i++)
                    for (j = 0; j < k; j++)
                        acb_set(acb_mat_entry(J, i, j), acb_mat_entry(Y, v[i], j));

                for (i = 0; i < k; i++)
                    acb_add(acb_mat_entry(J, i, i), acb_mat_entry(J, i, i), lambda, prec);
            }

            /* The correction for the frozen rows corresponds
               to the eigenvalue. */
            if (k == 1)
            {
                /* Just one eigenvalue. */
                acb_get_mag(eps, acb_mat_entry(Y, v[0], 0));
            }
            else
            {
                /* Inclusion of eigenvalues of lambda_approx I_k + V^T Y. */
                arb_mat_t M;
                arb_mat_init(M, k, k);

                /* Extract rows of Y corresponding to the eigenvalue correction. */
                for (i = 0; i < k; i++)
                {
                    for (j = 0; j < k; j++)
                    {
                        acb_get_mag(eps, acb_mat_entry(Y, v[i], j));
                        arf_set_mag(arb_midref(arb_mat_entry(M, i, j)), eps);
                    }
                }

                arb_mat_nonnegative_eig_bound(eps, M, prec);
                arb_mat_clear(M);
            }

            /* Error bound for eigenvalues. */
            acb_add_error_mag(lambda, eps);

            acb_mat_get_mid(X, X_approx);
            /* Error bounds for eigenvectors. */
            /* Update the not frozen rows of the eigenvectors. */
            for (i = 0; i < n - k; i++)
            {
                for (j = 0; j < k; j++)
                    acb_add(acb_mat_entry(X, u[i], j),
                            acb_mat_entry(X, u[i], j),
                            acb_mat_entry(Y, u[i], j), prec);
            }

            goto cleanup;
        }
    }

    /* We failed to find an enclosure. */
    acb_indeterminate(lambda);
    acb_mat_indeterminate(X);
    if (J != NULL)
        acb_mat_indeterminate(J);

cleanup:
    acb_mat_clear(R);
    acb_mat_clear(T);
    acb_mat_clear(Y);
    acb_mat_clear(Y0);
    acb_mat_clear(Yeps);
    acb_mat_clear(UY);
    acb_mat_clear(VY);
    mag_clear(eps);
    flint_free(u);
    flint_free(v);
}
