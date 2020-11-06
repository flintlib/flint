/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

void
ca_mat_set_jordan_blocks(ca_mat_t mat, const ca_vec_t lambda, slong num_blocks, slong * block_lambda, slong * block_size, ca_ctx_t ctx)
{
    slong i, j, n, count;

    n = ca_mat_nrows(mat);

    if (n != ca_mat_ncols(mat))
    {
        flint_printf("ca_mat_set_jordan_blocks: matrix must be square\n");
        flint_abort();
    }

    count = 0;
    for (i = 0; i < num_blocks; i++)
        count += block_size[i];

    if (count != n)
    {
        flint_printf("ca_mat_set_jordan_blocks: sum of block sizes does not agree with size of output matrix\n");
        flint_abort();
    }

    ca_mat_zero(mat, ctx);

    count = 0;
    for (i = 0; i < num_blocks; i++)
    {
        for (j = 0; j < block_size[i]; j++)
        {
            ca_set(ca_mat_entry(mat, count + j, count + j), ca_vec_entry(lambda, block_lambda[i]), ctx);

            if (j < block_size[i] - 1)
                ca_one(ca_mat_entry(mat, count + j, count + j + 1), ctx);
        }

        count += block_size[i];
    }
}

int
ca_mat_jordan_blocks(ca_vec_t lambda, slong * num_blocks, slong * block_lambda, slong * block_size, const ca_mat_t A, ca_ctx_t ctx)
{
    slong i, j, k, n;
    slong * ranks, * diagram;
    slong ranks_len, rank;
    ulong * exp;
    int success;

    n = ca_mat_nrows(A);

    if (n != ca_mat_ncols(A))
    {
        flint_printf("ca_mat_jordan_blocks: matrix must be square\n");
        flint_abort();
    }

    exp = flint_malloc(sizeof(ulong) * n);
    ranks = flint_malloc(sizeof(slong) * (n + 1));
    diagram = flint_malloc(sizeof(slong) * n);

    success = ca_mat_eigenvalues(lambda, exp, A, ctx);

    if (success)
    {
        *num_blocks = 0;

        for (i = 0; success && i < ca_vec_length(lambda, ctx); i++)
        {
            if (exp[i] == 1)
            {
                block_lambda[*num_blocks] = i;
                block_size[*num_blocks] = 1;
                *num_blocks += 1;
            }
            else
            {
                ca_mat_t B, C;

                ca_mat_init(B, n, n, ctx);
                ca_mat_init(C, n, n, ctx);

                for (j = 0; j < n; j++)
                    for (k = 0; k < n; k++)
                        if (j == k)
                            ca_sub(ca_mat_entry(B, j, j), ca_mat_entry(A, j, j), ca_vec_entry(lambda, i), ctx);
                        else
                            ca_set(ca_mat_entry(B, j, k), ca_mat_entry(A, j, k), ctx);

                ca_mat_set(C, B, ctx);

                success = ca_mat_rank(&rank, C, ctx);
                
                ranks_len = 2;
                ranks[0] = n;
                ranks[1] = rank;

                j = 0;
                while (success && (ranks[j] > ranks[j + 1] && ranks[j + 1] + exp[i] > n))
                {
                    ca_mat_mul(C, B, C, ctx);

                    success = ca_mat_rank(&rank, C, ctx);
                    ranks[ranks_len] = rank;
                    j++;
                    ranks_len++;
                }

                if (success)
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

                ca_mat_clear(B, ctx);
                ca_mat_clear(C, ctx);
            }
        }
    }

    flint_free(exp);
    flint_free(ranks);
    flint_free(diagram);

    return success;
}
