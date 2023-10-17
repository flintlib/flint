/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static void
sp2gz_reduce_delta(fmpz_mat_t res, fmpz_mat_t w, const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t delta, col, u;
    fmpz_t d;

    fmpz_mat_init(u, g, g);
    fmpz_mat_init(col, g, 1);
    fmpz_mat_window_init(delta, mat, 2 * g - 1, g, 2 * g, 2 * g);
    fmpz_init(d);

    fmpz_mat_transpose(col, delta);
    fmpz_mat_hnf_transform(col, u, col);
    fmpz_mat_inv(u, d, u);
    if (!fmpz_equal_si(d, 1))
    {
        flint_printf("(sp2gz_decompose) Error: not invertible\n");
        flint_abort();
    }
    sp2gz_block_diag(w, u);
    fmpz_mat_mul(res, w, mat);

    fmpz_mat_clear(u);
    fmpz_mat_clear(col);
    fmpz_mat_window_clear(delta);
    fmpz_clear(d);
}

static slong
sp2gz_reduce_gamma(fmpz_mat_t next, fmpz_mat_struct* vec, const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t cur, u;
    fmpz_t d, r;
    slong k;
    slong res = 0;
    int test = 0;

    /* Return 0 if gamma is already zero */
    for (k = 0; (k < g) && !test; k++)
    {
        if (!fmpz_is_zero(fmpz_mat_entry(mat, 2 * g - 1, k)))
        {
            test = 1;
        }
    }
    if (!test)
    {
        fmpz_mat_set(next, mat);
        return res;
    }

    fmpz_mat_init(cur, 2 * g, 2 * g);
    fmpz_mat_init(u, g, g);
    fmpz_init(d);
    fmpz_init(r);

    sp2gz_j(&vec[res]);
    fmpz_mat_mul(cur, mat, &vec[res]);
    res++;

    test = 0;
    for (k = 0; k < g - 1; k++)
    {
        if (!fmpz_is_zero(fmpz_mat_entry(cur, 2 * g - 1, k)))
        {
            test = 1;
        }
    }
    if (test)
    {
        sp2gz_reduce_delta(cur, &vec[res], cur);
        res++;
    }

    /* Set last row of u such that |last row of c + du| <= d/2 */
    for (k = 0; k < g; k++)
    {
        fmpz_set(d, fmpz_mat_entry(cur, 2 * g - 1, 2 * g - 1));
        fmpz_smod(r, fmpz_mat_entry(cur, 2 * g - 1, k), d);
        fmpz_sub(r, r, fmpz_mat_entry(cur, 2 * g - 1, k));
        fmpz_divexact(fmpz_mat_entry(u, g - 1, k), r, d);
        fmpz_set(fmpz_mat_entry(u, k, g - 1), fmpz_mat_entry(u, g - 1, k));
    }
    if (!fmpz_mat_is_zero(u))
    {
        sp2gz_j(&vec[res]);
        fmpz_mat_mul(cur, cur, &vec[res]);
        res++;

        fmpz_mat_neg(u, u);
        sp2gz_trig(&vec[res], u);
        fmpz_mat_mul(cur, cur, &vec[res]);
        res++;

        sp2gz_j(&vec[res]);
        fmpz_mat_mul(cur, cur, &vec[res]);
        res++;
    }

    flint_printf("start and new matrix:\n");
    fmpz_mat_print_pretty(mat);
    flint_printf("\n");
    fmpz_mat_print_pretty(cur);
    flint_printf("\n\n");
    fmpz_mat_set(next, cur);

    fmpz_mat_clear(cur);
    fmpz_mat_clear(u);
    fmpz_clear(d);
    fmpz_clear(r);
    return res;
}

static void
sp2gz_embed(fmpz_mat_t res, const fmpz_mat_t mat)
{
    slong j, k;
    slong g = sp2gz_dim(res);
    slong g1 = sp2gz_dim(mat);

    fmpz_mat_one(res);
    for (j = g - g1; j < g + g1; j++)
    {
        for (k = g - g1; k < g + g1; k++)
        {
            fmpz_set(fmpz_mat_entry(res, j, k), fmpz_mat_entry(mat, j, k));
        }
    }
}

fmpz_mat_struct* sp2gz_decompose(slong* nb, const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t cur;
    fmpz_mat_struct* gamma;
    fmpz_mat_struct* rec = NULL;
    fmpz_mat_struct* res;
    fmpz_mat_t w;
    slong nb_max = 0;
    slong nb_rec = 0;
    slong k, nb_gamma, add;

    fmpz_mat_init(cur, 2 * g, 2 * g);

    /* We need at most 5 * bits matrices to reduce gamma to zero */
    for (k = 0; k < g; k++)
    {
        nb_max = FLINT_MAX(nb_max, fmpz_bits(fmpz_mat_entry(mat, 2 * g - 1, k)));
        nb_max = 5 * nb_max + 1; /* for last delta reduction */
    }

    gamma = flint_malloc(nb_max * sizeof(fmpz_mat_struct));
    for (k = 0; k < nb_max; k++)
    {
        fmpz_mat_init(&gamma[k], 2 * g, 2 * g);
    }

    nb_gamma = 0;
    add = 1;
    fmpz_mat_set(cur, mat);
    while (add > 0)
    {
        add = sp2gz_reduce_gamma(cur, gamma + nb_gamma, cur);
        nb_gamma += add;
    }

    /* Reduce delta one last time before recursive call */
    sp2gz_reduce_delta(cur, &gamma[nb_gamma], cur);
    if (!fmpz_mat_is_one(&gamma[nb_gamma]))
    {
        nb_gamma++;
    }
    if (g > 1)
    {
        fmpz_mat_window_init(w, cur, 1, 1, g - 1, g - 1);
        rec = sp2gz_decompose(&nb_rec, w);
        fmpz_mat_window_clear(w);
    }

    /* Make final vector */
    *nb = nb_gamma + nb_rec;
    res = flint_malloc(*nb * sizeof(fmpz_mat_struct));
    for (k = 0; k < *nb; k++)
    {
        fmpz_mat_init(&res[k], 2 * g, 2 * g);
    }

    for (k = 0; k < nb_gamma; k++)
    {
        fmpz_mat_set(&res[k], &gamma[k]);
    }
    for (k = 0; k < nb_rec; k++)
    {
        sp2gz_embed(&res[nb_gamma + k], &rec[k]);
    }

    fmpz_mat_clear(cur);
    for (k = 0; k < nb_max; k++)
    {
        fmpz_mat_clear(&gamma[k]);
    }
    flint_free(gamma);
    for (k = 0; k < nb_rec; k++)
    {
        fmpz_mat_clear(&rec[k]);
    }
    flint_free(rec);
    return res;
}
