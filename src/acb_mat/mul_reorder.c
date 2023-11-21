/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

static void
copy_re_shallow(arb_mat_t X, const acb_mat_t A)
{
    slong M, N, i, j;
    M = arb_mat_nrows(X);
    N = arb_mat_ncols(X);
    for (i = 0; i < M; i++)
        for (j = 0; j < N; j++)
            *arb_mat_entry(X, i, j) = *acb_realref(acb_mat_entry(A, i, j));
}

static void
copy_im_shallow(arb_mat_t X, const acb_mat_t A)
{
    slong M, N, i, j;
    M = arb_mat_nrows(X);
    N = arb_mat_ncols(X);
    for (i = 0; i < M; i++)
        for (j = 0; j < N; j++)
            *arb_mat_entry(X, i, j) = *acb_imagref(acb_mat_entry(A, i, j));
}

static void
clear_shallow(arb_mat_t X)
{
    slong M, N, i, j;
    M = arb_mat_nrows(X);
    N = arb_mat_ncols(X);
    for (i = 0; i < M; i++)
        for (j = 0; j < N; j++)
            arb_init(arb_mat_entry(X, i, j));
}

/* todo: squaring optimizations */
void
acb_mat_mul_reorder(acb_mat_t C, const acb_mat_t A, const acb_mat_t B, slong prec)
{
    arb_mat_t X, Y, Z, W;
    slong M, N, P;
    slong i, j;

    M = acb_mat_nrows(A);
    N = acb_mat_ncols(A);
    P = acb_mat_ncols(B);

    if (acb_mat_is_real(A))
    {
        if (acb_mat_is_real(B))
        {
            arb_mat_init(X, M, N);
            arb_mat_init(Y, N, P);
            arb_mat_init(Z, M, P);

            copy_re_shallow(X, A);
            copy_re_shallow(Y, B);

            for (i = 0; i < M; i++)
                for (j = 0; j < P; j++)
                    arb_zero(acb_imagref(acb_mat_entry(C, i, j)));

            if (A == C || B == C)
            {
                arb_mat_mul(Z, X, Y, prec);

                for (i = 0; i < M; i++)
                    for (j = 0; j < P; j++)
                        arb_swap(acb_realref(acb_mat_entry(C, i, j)), acb_mat_entry(Z, i, j));
            }
            else
            {
                copy_re_shallow(Z, C);
                arb_mat_mul(Z, X, Y, prec);

                for (i = 0; i < M; i++)
                    for (j = 0; j < P; j++)
                        *acb_realref(acb_mat_entry(C, i, j)) = *arb_mat_entry(Z, i, j);

                clear_shallow(Z);
            }

            clear_shallow(X);
            clear_shallow(Y);

            arb_mat_clear(X);
            arb_mat_clear(Y);
            arb_mat_clear(Z);
        }
        else
        {
            arb_mat_init(X, M, N);
            arb_mat_init(Y, N, P);
            arb_mat_init(Z, M, P);

            /* (reA*reB, reA*imB) */

            copy_re_shallow(X, A);
            copy_re_shallow(Y, B);

            if (A == C || B == C)
            {
                arb_mat_t T;
                arb_mat_init(T, M, P);

                arb_mat_mul(T, X, Y, prec);
                copy_im_shallow(Y, B);
                arb_mat_mul(Z, X, Y, prec);

                for (i = 0; i < M; i++)
                    for (j = 0; j < P; j++)
                        arb_swap(acb_realref(acb_mat_entry(C, i, j)), acb_mat_entry(T, i, j));

                for (i = 0; i < M; i++)
                    for (j = 0; j < P; j++)
                        arb_swap(acb_imagref(acb_mat_entry(C, i, j)), acb_mat_entry(Z, i, j));

                arb_mat_clear(T);
            }
            else
            {
                copy_re_shallow(Z, C);
                arb_mat_mul(Z, X, Y, prec);

                for (i = 0; i < M; i++)
                    for (j = 0; j < P; j++)
                        *acb_realref(acb_mat_entry(C, i, j)) = *arb_mat_entry(Z, i, j);

                copy_im_shallow(Z, C);
                copy_im_shallow(Y, B);
                arb_mat_mul(Z, X, Y, prec);

                for (i = 0; i < M; i++)
                    for (j = 0; j < P; j++)
                        *acb_imagref(acb_mat_entry(C, i, j)) = *arb_mat_entry(Z, i, j);

                clear_shallow(Z);
            }

            clear_shallow(X);
            clear_shallow(Y);

            arb_mat_clear(X);
            arb_mat_clear(Y);
            arb_mat_clear(Z);
        }
    }
    else if (acb_mat_is_real(B))
    {
        arb_mat_init(X, M, N);
        arb_mat_init(Y, N, P);
        arb_mat_init(Z, M, P);

        /* (reA*reB, imA*reB) */

        copy_re_shallow(X, A);
        copy_re_shallow(Y, B);

        if (A == C || B == C)
        {
            arb_mat_t T;
            arb_mat_init(T, M, P);

            arb_mat_mul(T, X, Y, prec);
            copy_im_shallow(X, A);
            arb_mat_mul(Z, X, Y, prec);

            for (i = 0; i < M; i++)
                for (j = 0; j < P; j++)
                    arb_swap(acb_realref(acb_mat_entry(C, i, j)), acb_mat_entry(T, i, j));

            for (i = 0; i < M; i++)
                for (j = 0; j < P; j++)
                    arb_swap(acb_imagref(acb_mat_entry(C, i, j)), acb_mat_entry(Z, i, j));

            arb_mat_clear(T);
        }
        else
        {
            copy_re_shallow(Z, C);
            arb_mat_mul(Z, X, Y, prec);

            for (i = 0; i < M; i++)
                for (j = 0; j < P; j++)
                    *acb_realref(acb_mat_entry(C, i, j)) = *arb_mat_entry(Z, i, j);

            copy_im_shallow(Z, C);
            copy_im_shallow(X, A);
            arb_mat_mul(Z, X, Y, prec);

            for (i = 0; i < M; i++)
                for (j = 0; j < P; j++)
                    *acb_imagref(acb_mat_entry(C, i, j)) = *arb_mat_entry(Z, i, j);

            clear_shallow(Z);
        }

        clear_shallow(X);
        clear_shallow(Y);

        arb_mat_clear(X);
        arb_mat_clear(Y);
        arb_mat_clear(Z);
    }
    else
    {
        arb_mat_init(X, M, N);
        arb_mat_init(Y, N, P);
        arb_mat_init(Z, M, P);
        arb_mat_init(W, M, P);

        copy_re_shallow(X, A);
        copy_re_shallow(Y, B);
        arb_mat_mul(Z, X, Y, prec);

        copy_im_shallow(X, A);
        copy_im_shallow(Y, B);
        arb_mat_mul(W, X, Y, prec);

        if (A == C || B == C)
        {
            arb_mat_t T;
            arb_mat_init(T, M, P);
            arb_mat_sub(T, Z, W, prec);

            copy_re_shallow(X, A);
            arb_mat_mul(Z, X, Y, prec);

            copy_im_shallow(X, A);
            copy_re_shallow(Y, B);
            arb_mat_mul(W, X, Y, prec);

            for (i = 0; i < M; i++)
                for (j = 0; j < P; j++)
                    arb_swap(acb_realref(acb_mat_entry(C, i, j)),
                        arb_mat_entry(T, i, j));

            arb_mat_clear(T);

            for (i = 0; i < M; i++)
                for (j = 0; j < P; j++)
                    arb_add(acb_imagref(acb_mat_entry(C, i, j)),
                        arb_mat_entry(Z, i, j), arb_mat_entry(W, i, j), prec);
        }
        else
        {
            for (i = 0; i < M; i++)
                for (j = 0; j < P; j++)
                    arb_sub(acb_realref(acb_mat_entry(C, i, j)),
                        arb_mat_entry(Z, i, j), arb_mat_entry(W, i, j), prec);

            copy_re_shallow(X, A);
            arb_mat_mul(Z, X, Y, prec);

            copy_im_shallow(X, A);
            copy_re_shallow(Y, B);
            arb_mat_mul(W, X, Y, prec);

            for (i = 0; i < M; i++)
                for (j = 0; j < P; j++)
                    arb_add(acb_imagref(acb_mat_entry(C, i, j)),
                        arb_mat_entry(Z, i, j), arb_mat_entry(W, i, j), prec);
        }

        clear_shallow(X);
        clear_shallow(Y);

        arb_mat_clear(X);
        arb_mat_clear(Y);
        arb_mat_clear(Z);
        arb_mat_clear(W);
    }
}

