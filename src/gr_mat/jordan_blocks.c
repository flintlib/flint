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

int
gr_mat_set_jordan_blocks(gr_mat_t mat, const gr_vec_t lambda, slong num_blocks, slong * block_lambda, slong * block_size, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i, j, n, count;

    n = gr_mat_nrows(mat, ctx);

    if (n != gr_mat_ncols(mat, ctx))
    {
        /* matrix must be square */
        return GR_DOMAIN;
    }

    count = 0;
    for (i = 0; i < num_blocks; i++)
        count += block_size[i];

    if (count != n)
    {
        /* sum of block sizes does not agree with size of output matrix */
        return GR_DOMAIN;
    }

    status |= gr_mat_zero(mat, ctx);

    count = 0;
    for (i = 0; i < num_blocks; i++)
    {
        for (j = 0; j < block_size[i]; j++)
        {
            status |= gr_set(gr_mat_entry_ptr(mat, count + j, count + j, ctx), gr_vec_entry_srcptr(lambda, block_lambda[i], ctx), ctx);

            if (j < block_size[i] - 1)
                status |= gr_one(gr_mat_entry_ptr(mat, count + j, count + j + 1, ctx), ctx);
        }

        count += block_size[i];
    }

    return status;
}

int
gr_mat_jordan_blocks(gr_vec_t lambda, slong * num_blocks, slong * block_lambda, slong * block_size, const gr_mat_t A, gr_ctx_t ctx)
{
    slong i, j, k, n;
    slong * ranks, * diagram;
    slong ranks_len, rank;
    fmpz * exp;
    slong exp_i;
    int status = GR_SUCCESS;
    gr_ctx_t ZZ;
    gr_vec_t mult;

    n = gr_mat_nrows(A, ctx);

    if (n != gr_mat_ncols(A, ctx))
    {
        /* matrix must be square */
        return GR_DOMAIN;
    }

    ranks = flint_malloc(sizeof(slong) * (n + 1));
    diagram = flint_malloc(sizeof(slong) * n);

    gr_ctx_init_fmpz(ZZ);
    gr_vec_init(mult, 0, ZZ);

    status |= gr_mat_eigenvalues(lambda, mult, A, 0, ctx);

    if (status == GR_SUCCESS)
    {
        exp = mult->entries;

        *num_blocks = 0;

        for (i = 0; (status == GR_SUCCESS) && i < gr_vec_length(lambda, ctx); i++)
        {
            exp_i = exp[i];

            if (exp_i == 1)
            {
                block_lambda[*num_blocks] = i;
                block_size[*num_blocks] = 1;
                *num_blocks += 1;
            }
            else
            {
                gr_mat_t B, C;

                gr_mat_init(B, n, n, ctx);
                gr_mat_init(C, n, n, ctx);

                for (j = 0; j < n; j++)
                    for (k = 0; k < n; k++)
                        if (j == k)
                            status |= gr_sub(gr_mat_entry_ptr(B, j, j, ctx), gr_mat_entry_srcptr(A, j, j, ctx), gr_vec_entry_srcptr(lambda, i, ctx), ctx);
                        else
                            status |= gr_set(gr_mat_entry_ptr(B, j, k, ctx), gr_mat_entry_srcptr(A, j, k, ctx), ctx);

                status |= gr_mat_set(C, B, ctx);

                status |= gr_mat_rank(&rank, C, ctx);

                ranks_len = 2;
                ranks[0] = n;
                ranks[1] = rank;

                j = 0;
                while ((status == GR_SUCCESS) && (ranks[j] > ranks[j + 1] && ranks[j + 1] + exp_i > n))
                {
                    status |= gr_mat_mul(C, B, C, ctx);
                    status |= gr_mat_rank(&rank, C, ctx);
                    ranks[ranks_len] = rank;
                    j++;
                    ranks_len++;
                }

                if (status == GR_SUCCESS)
                {
                    /* Ferrer's diagram of an integer partition */
                    for (j = 0; j < ranks_len - 1; j++)
                        diagram[j] = ranks[j] - ranks[j + 1];

                    /* Transpose Ferrer's diagram */
                    for (j = 1; j <= diagram[0]; j++)
                    {
                        slong c = 0;

                        for (k = 0; k < ranks_len - 1; k++)
                            c += (diagram[k] >= j);

                        block_lambda[*num_blocks] = i;
                        block_size[*num_blocks] = c;
                        *num_blocks += 1;
                    }
                }

                gr_mat_clear(B, ctx);
                gr_mat_clear(C, ctx);
            }
        }
    }

    gr_vec_clear(mult, ZZ);
    gr_ctx_clear(ZZ);

    flint_free(ranks);
    flint_free(diagram);

    return status;
}
