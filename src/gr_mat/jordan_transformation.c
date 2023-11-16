/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_mat.h"

/* The algorithm is taken from jordan_form() in Sage */

/* Find a row vector v in the row span of V but not in the row span of W.
   Returns the index of such a vector.
   Returns -1 on failure. */
static slong
vector_in_difference(const gr_mat_t V, const gr_mat_t W, gr_ctx_t ctx)
{
    gr_mat_t U;
    gr_ptr v;
    gr_ptr t, u;
    slong i, j, k, l, n, found, rank;
    truth_t is_zero;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (gr_mat_nrows(V, ctx) == 0)
        return -1;

    if (gr_mat_nrows(W, ctx) == 0)
        return 0;

    n = gr_mat_ncols(W, ctx);
    found = -1;

    gr_mat_init(U, gr_mat_nrows(W, ctx), n, ctx);
    GR_TMP_INIT_VEC(v, n, ctx);
    GR_TMP_INIT2(t, u, ctx);

    if (gr_mat_rref(&rank, U, W, ctx) == GR_SUCCESS)
    {
        for (i = 0; i < gr_mat_nrows(V, ctx); i++)  /* for candidate v in V */
        {
            /* Copy of v for reduction */
            GR_MUST_SUCCEED(_gr_vec_set(v, gr_mat_entry_srcptr(V, i, 0, ctx), n, ctx));

            /* Reduce v by the rows in W */
            for (j = 0; j < rank; j++)         /* for each w in W */
            {
                for (k = 0; k < n; k++)  /* find pivot element in w */
                {
                    is_zero = gr_is_zero(gr_mat_entry_srcptr(U, j, k, ctx), ctx);

                    if (is_zero == T_UNKNOWN)
                        goto cleanup;

                    /* reduce by this row */
                    if (is_zero == T_FALSE)
                    {
                        status |= gr_div(t, GR_ENTRY(v, k, sz), gr_mat_entry_srcptr(U, j, k, ctx), ctx);

                        for (l = 0; l < n; l++)
                        {
                            if (l == k)
                                status |= gr_zero(GR_ENTRY(v, l, sz), ctx);
                            else
                            {
                                status |= gr_mul(u, t, gr_mat_entry_srcptr(U, j, l, ctx), ctx);
                                status |= gr_sub(GR_ENTRY(v, l, sz), GR_ENTRY(v, l, sz), u, ctx);
                            }
                        }

                        break;
                    }
                }
            }

            if (status != GR_SUCCESS)
                goto cleanup;

            is_zero = _gr_vec_is_zero(v, n, ctx);

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
    gr_mat_clear(U, ctx);
    GR_TMP_CLEAR_VEC(v, n, ctx);
    GR_TMP_CLEAR2(t, u, ctx);

    return found;
}

/* note: doesn't support aliasing */
static int
_gr_mat_mul_vec(gr_ptr res, const gr_mat_t mat, gr_srcptr v, gr_ctx_t ctx)
{
    slong i;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    for (i = 0; i < gr_mat_nrows(mat, ctx); i++)
    {
        status |= _gr_vec_dot(GR_ENTRY(res, i, sz), NULL, 0, gr_mat_entry_srcptr(mat, i, 0, ctx), v, gr_mat_ncols(mat, ctx), ctx);
    }

    return status;
}

static void
gr_mat_transpose_resize(gr_mat_t B, const gr_mat_t A, gr_ctx_t ctx)
{
    gr_mat_t T;
    gr_mat_init(T, gr_mat_ncols(A, ctx), gr_mat_nrows(A, ctx), ctx);
    GR_MUST_SUCCEED(gr_mat_transpose(T, A, ctx));
    gr_mat_swap(B, T, ctx);
    gr_mat_clear(T, ctx);
}

static int
gr_mat_pow_ui_binexp(gr_mat_t res, const gr_mat_t A, ulong n, gr_ctx_t ctx)
{
    gr_ctx_t mctx;
    int status;

    gr_ctx_init_matrix_ring(mctx, ctx, gr_mat_nrows(A, ctx));
    status = gr_pow_ui(res, A, n, mctx);
    gr_ctx_clear(mctx);

    return status;
}

/* todo: special-case num_lambda == n (diagonalization) */
/* todo: also, for each vector */
/* todo: aliasing */
int
gr_mat_jordan_transformation(gr_mat_t mat, const gr_vec_t lambda, slong num_blocks, slong * block_lambda, slong * block_size, const gr_mat_t A, gr_ctx_t ctx)
{
    slong num_lambda, i, j, k, l, m, n, output_block, column_offset;
    slong *sizes, *counts;
    gr_mat_t B, Y, V1, V2, V1ker, V2ker;
    slong size, num_sizes, count, y_rows, v_index;
    int * written;
    int status;

    n = gr_mat_nrows(A, ctx);

    if (n == 0)
        return GR_SUCCESS;

    num_lambda = gr_vec_length(lambda, ctx);

    /* Special-case diagonalization with distinct eigenvalues. */
    /* TODO: also special-case eigenvalues with multiplicity 1 below. */
    if (num_lambda == n)
    {
        status = GR_SUCCESS;

        gr_mat_init(B, n, n, ctx);
        gr_mat_init(Y, 0, 0, ctx);

        for (i = 0; i < n; i++)
        {
            /* B = A - lambda */
            for (j = 0; j < n; j++)
                for (k = 0; k < n; k++)
                    if (j == k)
                        status |= gr_sub(gr_mat_entry_ptr(B, j, j, ctx), gr_mat_entry_srcptr(A, j, j, ctx), gr_vec_entry_srcptr(lambda, block_lambda[i], ctx), ctx);
                    else
                        status |= gr_set(gr_mat_entry_ptr(B, j, k, ctx), gr_mat_entry_srcptr(A, j, k, ctx), ctx);

            status |= gr_mat_nullspace(Y, B, ctx);

            if (status != GR_SUCCESS)
                break;

            if (gr_mat_ncols(Y, ctx) != 1)
            {
                /* should not happen */
                status = GR_UNABLE;
            }

            for (j = 0; j < n; j++)
                status |= gr_set(gr_mat_entry_ptr(mat, j, i, ctx), gr_mat_entry_srcptr(Y, j, 0, ctx), ctx);
        }

        gr_mat_clear(B, ctx);
        gr_mat_clear(Y, ctx);

        return status;
    }

    sizes = flint_malloc(sizeof(slong) * n);
    counts = flint_malloc(sizeof(slong) * n);
    written = flint_calloc(num_blocks, sizeof(int));
    gr_mat_init(B, n, n, ctx);
    gr_mat_init(Y, 0, n, ctx);
    gr_mat_init(V1, n, n, ctx);
    gr_mat_init(V2, n, n, ctx);
    gr_mat_init(V1ker, 0, 0, ctx);
    gr_mat_init(V2ker, 0, 0, ctx);

    status = GR_SUCCESS;

    for (i = 0; i < num_lambda; i++)
    {
        /* B = A - lambda */
        for (j = 0; j < n; j++)
            for (k = 0; k < n; k++)
                if (j == k)
                    status |= gr_sub(gr_mat_entry_ptr(B, j, j, ctx), gr_mat_entry_srcptr(A, j, j, ctx), gr_vec_entry_srcptr(lambda, i, ctx), ctx);
                else
                    status |= gr_set(gr_mat_entry_ptr(B, j, k, ctx), gr_mat_entry_srcptr(A, j, k, ctx), ctx);

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
        gr_mat_clear(Y, ctx);
        gr_mat_init(Y, n, n, ctx);
        y_rows = 0;

        /* Iterate over (block size, count) */
        for (j = 0; j < num_sizes; j++)
        {
            size = sizes[j];
            count = counts[j];

            /* should not happen */
            if (size == 0)
            {
                status = GR_UNABLE;
                goto cleanup;
            }

            /* Find elements in ker(B^size) - ker(B^(size-1)) */
            status |= gr_mat_pow_ui_binexp(V2, B, size - 1, ctx);
            status |= gr_mat_mul(V1, B, V2, ctx);

            status |= gr_mat_nullspace(V1ker, V1, ctx);
            if (status != GR_SUCCESS)
                goto cleanup;

            status |= gr_mat_nullspace(V2ker, V2, ctx);
            if (status != GR_SUCCESS)
                goto cleanup;

            /* rows instead of columns */
            gr_mat_transpose_resize(V1ker, V1ker, ctx);
            gr_mat_transpose_resize(V2ker, V2ker, ctx);

            for (k = 0; k < count; k++)
            {
                /* Concatenate V2ker with Y */
                gr_mat_t V2kerY;

                gr_mat_init(V2kerY, gr_mat_nrows(V2ker, ctx) + y_rows, n, ctx);
                for (m = 0; m < gr_mat_nrows(V2ker, ctx); m++)
                    status |= _gr_vec_set(gr_mat_entry_ptr(V2kerY, m, 0, ctx), gr_mat_entry_srcptr(V2ker, m, 0, ctx), n, ctx);
                for (m = 0; m < y_rows; m++)
                    status |= _gr_vec_set(gr_mat_entry_ptr(V2kerY, gr_mat_nrows(V2ker, ctx) + m, 0, ctx), gr_mat_entry_srcptr(Y, m, 0, ctx), n, ctx);

                v_index = vector_in_difference(V1ker, V2kerY, ctx);

                gr_mat_clear(V2kerY, ctx);

                if (v_index == -1)
                {
                    status = GR_UNABLE;
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
                status |= _gr_vec_set(gr_mat_entry_ptr(Y, y_rows, 0, ctx), gr_mat_entry_srcptr(V1ker, v_index, 0, ctx), n, ctx);
                for (m = 1; m < size; m++)
                {
                    status |= _gr_mat_mul_vec(gr_mat_entry_ptr(Y, y_rows + m, 0, ctx), B, gr_mat_entry_srcptr(Y, y_rows + m - 1, 0, ctx), ctx);
                }
                y_rows += size;

                /* Insert chain in the output matrix */
                for (m = 0; m < size; m++)
                {
                    for (l = 0; l < n; l++)
                        status |= gr_set(gr_mat_entry_ptr(mat, l, column_offset + m, ctx), gr_mat_entry_srcptr(Y, y_rows - 1 - m, l, ctx), ctx);
                }
            }
        }
    }

cleanup:
    flint_free(sizes);
    flint_free(counts);
    gr_mat_clear(B, ctx);
    gr_mat_clear(Y, ctx);
    gr_mat_clear(V1, ctx);
    gr_mat_clear(V2, ctx);
    gr_mat_clear(V1ker, ctx);
    gr_mat_clear(V2ker, ctx);
    flint_free(written);

    return status;
}
