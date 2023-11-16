/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

truth_t
ca_mat_diagonalization_precomp(ca_mat_t D, ca_mat_t P, const ca_mat_t A, const ca_vec_t eigenvalues, const ulong * am, ca_ctx_t ctx)
{
    int success;
    truth_t result;
    ca_mat_t AIe, b;
    slong i, j, k, n;
    slong nullity, added;

    n = ca_mat_nrows(A);

    ca_mat_init(AIe, n, n, ctx);
    ca_mat_init(b, 0, 0, ctx);

    result = T_TRUE;

    ca_mat_zero(D, ctx);

    added = 0;
    for (i = 0; i < ca_vec_length(eigenvalues, ctx); i++)
    {
        ca_mat_set(AIe, A, ctx);
        for (j = 0; j < n; j++)
            ca_sub(ca_mat_entry(AIe, j, j), ca_mat_entry(AIe, j, j), ca_vec_entry(eigenvalues, i), ctx);

        success = ca_mat_right_kernel(b, AIe, ctx);
        if (!success)
        {
            result = T_UNKNOWN;
            break;
        }

        nullity = ca_mat_ncols(b);

        if (nullity != am[i])
        {
            result = T_FALSE;
            break;
        }

        for (j = 0; j < am[i]; j++)
        {
            ca_set(ca_mat_entry(D, added + j, added + j), ca_vec_entry(eigenvalues, i), ctx);

            for (k = 0; k < n; k++)
                ca_set(ca_mat_entry(P, k, added + j), ca_mat_entry(b, k, j), ctx);
        }

        added += am[i];
    }

    ca_mat_clear(AIe, ctx);
    ca_mat_clear(b, ctx);

    return result;
}

truth_t
ca_mat_diagonalization(ca_mat_t D, ca_mat_t P, const ca_mat_t A, ca_ctx_t ctx)
{
    truth_t result;
    ca_vec_t eigenvalues;
    ulong * am;
    slong n;

    if (!ca_mat_is_square(A))
        return T_FALSE;

    n = ca_mat_nrows(A);

    am = flint_malloc(sizeof(ulong) * n);
    ca_vec_init(eigenvalues, 0, ctx);

    if (ca_mat_eigenvalues(eigenvalues, am, A, ctx))
    {
        result = ca_mat_diagonalization_precomp(D, P, A, eigenvalues, am, ctx);
    }
    else
    {
        result = T_UNKNOWN;
    }

    ca_vec_clear(eigenvalues, ctx);
    flint_free(am);

    return result;
}
