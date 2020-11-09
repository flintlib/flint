/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

truth_t
_ca_vec_check_is_zero(ca_srcptr vec, slong len, ca_ctx_t ctx)
{
    slong i;
    int have_unknown;
    truth_t is_zero;

    have_unknown = 0;
    for (i = 0; i < len; i++)
    {
        is_zero = ca_check_is_zero(vec + i, ctx);

        if (is_zero == T_FALSE)
            return T_FALSE;

        if (is_zero == T_UNKNOWN)
            have_unknown = 1;
    }

    if (have_unknown)
        return T_UNKNOWN;
    else
        return T_TRUE;
}

/* The algorithm is taken from jordan_form() in Sage */

/* Find a row vector v in the row span of V but not in the row span of W.
   Returns the index of such a vector.
   Returns -1 on failure. */
static slong
vector_in_difference(const ca_mat_t V, const ca_mat_t W, ca_ctx_t ctx)
{
    ca_mat_t U;
    ca_ptr v;
    ca_t t, u;
    slong i, j, k, l, n, found, rank;
    truth_t is_zero;

    if (ca_mat_nrows(V) == 0)
        return -1;

    if (ca_mat_nrows(W) == 0)
        return 0;

    n = ca_mat_ncols(W);
    found = -1;

    ca_mat_init(U, ca_mat_nrows(W), n, ctx);
    v = _ca_vec_init(n, ctx);
    ca_init(t, ctx);
    ca_init(u, ctx);

    if (ca_mat_rref(&rank, U, W, ctx))
    {
        for (i = 0; i < ca_mat_nrows(V); i++)  /* for candidate v in V */
        {
            /* Copy of v for reduction */
            _ca_vec_set(v, ca_mat_entry(V, i, 0), n, ctx);

            /* Reduce v by the rows in W */
            for (j = 0; j < rank; j++)         /* for each w in W */
            {
                for (k = 0; k < n; k++)  /* find pivot element in w */
                {
                    is_zero = ca_check_is_zero(ca_mat_entry(U, j, k), ctx);

                    if (is_zero == T_UNKNOWN)
                        goto cleanup;

                    /* reduce by this row */
                    if (is_zero == T_FALSE)
                    {
                        ca_div(t, v + k, ca_mat_entry(U, j, k), ctx);

                        for (l = 0; l < n; l++)
                        {
                            if (l == k)
                                ca_zero(v + l, ctx);
                            else
                            {
                                ca_mul(u, t, ca_mat_entry(U, j, l), ctx);
                                ca_sub(v + l, v + l, u, ctx);
                            }
                        }

                        break;
                    }
                }
            }

            is_zero = _ca_vec_check_is_zero(v, n, ctx);

            if (is_zero == T_UNKNOWN)
                goto cleanup;

            if (is_zero == T_FALSE)
            {
                found = i;
                break;
            }
        }
    }


cleanup:
    ca_mat_clear(U, ctx);
    _ca_vec_clear(v, n, ctx);
    ca_clear(t, ctx);
    ca_clear(u, ctx);

    return found;
}

/* note: doesn't support aliasing */
static void
_ca_mat_mul_vec(ca_ptr res, const ca_mat_t mat, ca_srcptr v, ca_ctx_t ctx)
{
    slong i;

    for (i = 0; i < ca_mat_nrows(mat); i++)
    {
        ca_dot(res + i, NULL, 0, ca_mat_entry(mat, i, 0), 1, v, 1, ca_mat_ncols(mat), ctx);
    }
}

static void
ca_mat_transpose_resize(ca_mat_t B, const ca_mat_t A, ca_ctx_t ctx)
{
    ca_mat_t T;
    ca_mat_init(T, ca_mat_ncols(A), ca_mat_nrows(A), ctx);
    ca_mat_transpose(T, A, ctx);
    ca_mat_swap(B, T, ctx);
    ca_mat_clear(T, ctx);
}

/* todo: special-case num_lambda == n (diagonalization) */
/* todo: also, for each vector */
/* todo: aliasing */
int
ca_mat_jordan_transformation(ca_mat_t mat, const ca_vec_t lambda, slong num_blocks, slong * block_lambda, slong * block_size, const ca_mat_t A, ca_ctx_t ctx)
{
    slong num_lambda, i, j, k, l, m, n, output_block, column_offset;
    slong *sizes, *counts;
    ca_mat_t B, Y, V1, V2, V1ker, V2ker;
    slong size, num_sizes, count, y_rows, v_index;
    int * written;
    int success;

    n = ca_mat_nrows(A);

    if (n == 0)
        return 1;

    num_lambda = ca_vec_length(lambda, ctx);

    /* Special-case diagonalization with distinct eigenvalues. */
    /* TODO: also special-case eigenvalues with multiplicity 1 below. */
    if (num_lambda == n)
    {
        success = 1;

        ca_mat_init(B, n, n, ctx);
        ca_mat_init(Y, 0, 0, ctx);

        for (i = 0; i < n; i++)
        {
            /* B = A - lambda */
            for (j = 0; j < n; j++)
                for (k = 0; k < n; k++)
                    if (j == k)
                        ca_sub(ca_mat_entry(B, j, j), ca_mat_entry(A, j, j), ca_vec_entry(lambda, block_lambda[i]), ctx);
                    else
                        ca_set(ca_mat_entry(B, j, k), ca_mat_entry(A, j, k), ctx);


            success = ca_mat_right_kernel(Y, B, ctx);


            if (!success)
                break;

            if (ca_mat_ncols(Y) != 1)
                abort();

            for (j = 0; j < n; j++)
                ca_set(ca_mat_entry(mat, j, i), ca_mat_entry(Y, j, 0), ctx);
        }

        ca_mat_clear(B, ctx);
        ca_mat_clear(Y, ctx);

        return success;
    }

    sizes = flint_malloc(sizeof(slong) * n);
    counts = flint_malloc(sizeof(slong) * n);
    written = flint_calloc(num_blocks, sizeof(int));
    ca_mat_init(B, n, n, ctx);
    ca_mat_init(Y, 0, n, ctx);
    ca_mat_init(V1, n, n, ctx);
    ca_mat_init(V2, n, n, ctx);
    ca_mat_init(V1ker, 0, 0, ctx);
    ca_mat_init(V2ker, 0, 0, ctx);

    success = 1;

    for (i = 0; i < num_lambda; i++)
    {
        /* B = A - lambda */
        for (j = 0; j < n; j++)
            for (k = 0; k < n; k++)
                if (j == k)
                    ca_sub(ca_mat_entry(B, j, j), ca_mat_entry(A, j, j), ca_vec_entry(lambda, i), ctx);
                else
                    ca_set(ca_mat_entry(B, j, k), ca_mat_entry(A, j, k), ctx);

        /* Group blocks as (block size, count) */
        /* Todo: does this need to be sorted? */
        num_sizes = 0;
        for (j = 0; j < num_blocks; j++)
        {
            if (block_lambda[j] == i)
            {
                for (k = 0; k < num_sizes; k++)
                {
                    if (sizes[k] == block_size[j])
                    {
                        counts[k]++;
                        break;
                    }
                }

                if (k == num_sizes)
                {
                    sizes[num_sizes] = block_size[j];
                    counts[num_sizes] = 1;
                    num_sizes++;
                }
            }
        }

        /* Y = matrix of row vectors spanning the Jordan chains for
           this eigenvalue. */
        ca_mat_clear(Y, ctx);
        ca_mat_init(Y, n, n, ctx);
        y_rows = 0;

        /* Iterate over (block size, count) */
        for (j = 0; j < num_sizes; j++)
        {
            size = sizes[j];
            count = counts[j];

            if (size == 0)
                flint_abort();

            /* Find elements in ker(B^size) - ker(B^(size-1)) */
            ca_mat_pow_ui_binexp(V2, B, size - 1, ctx);
            ca_mat_mul(V1, B, V2, ctx);

            success = ca_mat_right_kernel(V1ker, V1, ctx);
            if (!success)
                goto cleanup;

            success = ca_mat_right_kernel(V2ker, V2, ctx);
            if (!success)
                goto cleanup;

            /* rows instead of columns */
            ca_mat_transpose_resize(V1ker, V1ker, ctx);
            ca_mat_transpose_resize(V2ker, V2ker, ctx);

            for (k = 0; k < count; k++)
            {
                /* Concatenate V2ker with Y */
                ca_mat_t V2kerY;

                ca_mat_init(V2kerY, ca_mat_nrows(V2ker) + y_rows, n, ctx);
                for (m = 0; m < ca_mat_nrows(V2ker); m++)
                    _ca_vec_set(ca_mat_entry(V2kerY, m, 0), ca_mat_entry(V2ker, m, 0), n, ctx);
                for (m = 0; m < y_rows; m++)
                    _ca_vec_set(ca_mat_entry(V2kerY, ca_mat_nrows(V2ker) + m, 0), ca_mat_entry(Y, m, 0), n, ctx);

                v_index = vector_in_difference(V1ker, V2kerY, ctx);

                ca_mat_clear(V2kerY, ctx);

                if (v_index == -1)
                {
                    success = 0;
                    goto cleanup;
                }

                /* Position in output matrix to write chain. */
                column_offset = 0;
                output_block = 0;
                while (1)
                {
                    if (block_lambda[output_block] == i && block_size[output_block] == size && !written[output_block])
                    {
                        written[output_block] = 1;
                        break;
                    }
                    else
                    {
                        column_offset += block_size[output_block];
                        output_block++;
                    }
                }

                /* chain = [v, B*v, ..., B^(size-1)*v] */
                _ca_vec_set(ca_mat_entry(Y, y_rows, 0), ca_mat_entry(V1ker, v_index, 0), n, ctx);
                for (m = 1; m < size; m++)
                {
                    _ca_mat_mul_vec(ca_mat_entry(Y, y_rows + m, 0), B, ca_mat_entry(Y, y_rows + m - 1, 0), ctx);
                }
                y_rows += size;

                /* Insert chain in the output matrix */
                for (m = 0; m < size; m++)
                {
                    for (l = 0; l < n; l++)
                        ca_set(ca_mat_entry(mat, l, column_offset + m), ca_mat_entry(Y, y_rows - 1 - m, l), ctx);
                }
            }
        }
    }

cleanup:
    flint_free(sizes);
    flint_free(counts);
    ca_mat_clear(B, ctx);
    ca_mat_clear(Y, ctx);
    ca_mat_clear(V1, ctx);
    ca_mat_clear(V2, ctx);
    ca_mat_clear(V1ker, ctx);
    ca_mat_clear(V2ker, ctx);
    flint_free(written);

    return success;
}
