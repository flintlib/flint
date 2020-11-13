/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

/* todo: use dot products */
int
_ca_mat_charpoly_danilevsky_inplace(ca_ptr p, ca_mat_t A, ca_ctx_t ctx)
{
    slong n, n_input;
    slong i, j, k;
    ca_ptr V, W, T;
    ca_ptr t, b;
    ca_t c, h;
    int success;
    truth_t is_zero;
    slong plen;

    n = n_input = ca_mat_nrows(A);

    if (n == 0)
    {
        ca_one(p, ctx);
        return 1;
    }

    if (n == 1)
    {
        ca_neg(p + 0, ca_mat_entry(A, 0, 0), ctx);
        ca_one(p + 1, ctx);
        return 1;
    }

    ca_init(c, ctx);

    success = 1;

    i = 1;
    ca_init(h, ctx);
    ca_one(p, ctx);
    plen = 1;

    t = _ca_vec_init(n + 1, ctx);
    b = _ca_vec_init(n + 1, ctx);
    V = _ca_vec_init(n, ctx);
    W = _ca_vec_init(n, ctx);
    T = _ca_vec_init(n, ctx);

    while (i < n)
    {
        ca_set(h, ca_mat_entry(A, n - i, n - i - 1), ctx);

        while (1)
        {
            is_zero = ca_check_is_zero(h, ctx);
            if (is_zero == T_FALSE)
                break;
            if (is_zero == T_UNKNOWN)
            {
                success = 0;
                goto cleanup;
            }

            k = 1;
            while (k < n - i)
            {
                is_zero = ca_check_is_zero(ca_mat_entry(A, n - i, n - i - k - 1), ctx);
                if (is_zero == T_FALSE)
                    break;
                if (is_zero == T_UNKNOWN)
                {
                    success = 0;
                    goto cleanup;
                }
                k++;
            }

            if (k == n - i)
            {
                ca_one(b + i, ctx);
                for (k = 1; k <= i; k++)
                    ca_neg(b + k - 1, ca_mat_entry(A, n - i, n - k), ctx);
                _ca_poly_mul(t, p, plen, b, i + 1, ctx);
                plen += i;
                _ca_vec_swap(p, t, plen, ctx);

                n -= i;
                i = 1;

                if (n == 1)
                {
                    ca_one(b + 1, ctx);
                    ca_neg(b, ca_mat_entry(A, 0, 0), ctx);
                    _ca_poly_mul(t, p, plen, b, 2, ctx);
                    plen += 1;
                    _ca_vec_swap(p, t, plen, ctx);
                    goto cleanup;
                }
            }
            else
            {
                ca_ptr ptr;

                ptr = A->rows[n - i - k - 1];
                A->rows[n - i - k - 1] = A->rows[n - i - 1];
                A->rows[n - i - 1] = ptr;

                for (j = 1; j <= n - i + 1; j++)
                {
                    ca_swap(ca_mat_entry(A, j - 1, n - i - k - 1),
                                  ca_mat_entry(A, j - 1, n - i - 1), ctx);
                }
            }

            ca_set(h, ca_mat_entry(A, n - i, n - i - 1), ctx);
        }
      
        ca_neg(h, h, ctx);
        ca_inv(h, h, ctx);
      
        for (j = 1; j <= n; j++)
        {
            ca_mul(V + j - 1, ca_mat_entry(A, n - i, j - 1), h, ctx);
            ca_set(W + j - 1, ca_mat_entry(A, n - i, j - 1), ctx);
        } 

        ca_neg(h, h, ctx);
      
        for (j = 1; j <= n - i; j++)
        {
            for (k = 1; k <= n - i - 1; k++)
            {
                ca_mul(c, ca_mat_entry(A, j - 1, n - i - 1), V + k - 1, ctx);
                ca_add(ca_mat_entry(A, j - 1, k - 1), 
                              ca_mat_entry(A, j - 1, k - 1), c, ctx);
            }

            for (k = n - i + 1; k <= n; k++)
            {
                ca_mul(c, ca_mat_entry(A, j - 1, n - i - 1), V + k - 1, ctx);
                ca_add(ca_mat_entry(A, j - 1, k - 1), 
                              ca_mat_entry(A, j - 1, k - 1), c, ctx);
            }

            ca_mul(ca_mat_entry(A, j - 1, n - i - 1),
                           ca_mat_entry(A, j - 1, n - i - 1), h, ctx);
        }

        for (j = 1; j <= n - i - 1; j++)
        {
            ca_mul(ca_mat_entry(A, n - i - 1, j - 1), 
                           ca_mat_entry(A, n - i - 1, j - 1), W + n - i - 1, ctx);

            for (k = 1; k < n - i; k++)
            {
                ca_mul(c, ca_mat_entry(A, k - 1, j - 1), W + k - 1, ctx);
                ca_add(ca_mat_entry(A, n - i - 1, j - 1),
                              ca_mat_entry(A, n - i - 1, j - 1), c, ctx);
            }
        }

        for (j = n - i; j <= n - 1; j++)
        {
            ca_mul(ca_mat_entry(A, n - i - 1, j - 1), 
                           ca_mat_entry(A, n - i - 1, j - 1), W + n - i - 1, ctx);

            for (k = 1; k < n - i; k++)
            {
                ca_mul(c, ca_mat_entry(A, k - 1, j - 1), W + k - 1, ctx);
                ca_add(ca_mat_entry(A, n - i - 1, j - 1),
                              ca_mat_entry(A, n - i - 1, j - 1), c, ctx);
            }

            ca_add(ca_mat_entry(A, n - i - 1, j - 1),
                              ca_mat_entry(A, n - i - 1, j - 1), W + j, ctx);
        }

        ca_mul(ca_mat_entry(A, n - i - 1, n - 1), 
                        ca_mat_entry(A, n - i - 1, n - 1), W + n - i - 1, ctx);

        for (k = 1; k < n - i; k++)
        {
            ca_mul(c, ca_mat_entry(A, k - 1, n - 1), W + k - 1, ctx);
            ca_add(ca_mat_entry(A, n - i - 1, n - 1),
                           ca_mat_entry(A, n - i - 1, n - 1), c, ctx);
        }

        i++;
    }

    ca_one(b + n, ctx);
    for (i = 1; i <= n; i++)
        ca_neg(b + i - 1, ca_mat_entry(A, 0, n - i), ctx);
    _ca_poly_mul(t, p, plen, b, n + 1, ctx);
    plen += n;
    _ca_vec_swap(p, t, plen, ctx);

cleanup:
  
    ca_clear(c, ctx);
    ca_clear(h, ctx);
    _ca_vec_clear(t, n_input + 1, ctx);
    _ca_vec_clear(b, n_input + 1, ctx);
    _ca_vec_clear(T, n_input, ctx);
    _ca_vec_clear(V, n_input, ctx);
    _ca_vec_clear(W, n_input, ctx);

    return success;
}

int
_ca_mat_charpoly_danilevsky(ca_ptr cp, const ca_mat_t A, ca_ctx_t ctx)
{
    ca_mat_t T;
    int success;

    ca_mat_init(T, ca_mat_nrows(A), ca_mat_nrows(A), ctx);
    ca_mat_set(T, A, ctx);
    success = _ca_mat_charpoly_danilevsky_inplace(cp, T, ctx);
    ca_mat_clear(T, ctx);
    return success;
}

int
ca_mat_charpoly_danilevsky(ca_poly_t cp, const ca_mat_t mat, ca_ctx_t ctx)
{
    ca_poly_fit_length(cp, mat->r + 1, ctx);
    _ca_poly_set_length(cp, mat->r + 1, ctx);
    return _ca_mat_charpoly_danilevsky(cp->coeffs, mat, ctx);
}
