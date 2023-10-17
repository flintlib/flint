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
    fmpz_mat_t delta, diag, u, v;
    slong k;

    fmpz_mat_init(u, g, g);
    fmpz_mat_init(v, g, g);
    fmpz_mat_init(diag, g, g);
    fmpz_mat_window_init(delta, mat, g, g, 2 * g, 2 * g);

    for (k = 0; k < g; k++)
    {
        fmpz_one(fmpz_mat_entry(diag, k, g - 1 - k));
    }
    fmpz_mat_transpose(v, delta);
    fmpz_mat_mul(v, v, diag);
    fmpz_mat_hnf_transform(v, u, v);
    fmpz_mat_mul(u, diag, u);

    sp2gz_block_diag(w, u);
    sp2gz_inv(w, w);
    fmpz_mat_mul(res, mat, w);

    fmpz_mat_clear(u);
    fmpz_mat_clear(v);
    fmpz_mat_clear(diag);
    fmpz_mat_window_clear(delta);
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

    /* Reduce last line of gamma */
    sp2gz_j(cur);
    fmpz_mat_mul(cur, mat, cur);
    sp2gz_reduce_delta(cur, &vec[res], cur);
    fmpz_mat_transpose(&vec[res], &vec[res]);
    sp2gz_inv(&vec[res], &vec[res]);
    fmpz_mat_mul(cur, mat, &vec[res]);
    if (!fmpz_mat_is_one(&vec[res]))
    {
        res++;
    }

    /* Reduce last line of delta mod d */
    fmpz_set(d, fmpz_mat_entry(cur, 2 * g - 1, g - 1));
    for (k = 0; k < g; k++)
    {
        fmpz_smod(r, fmpz_mat_entry(cur, 2 * g - 1, g + k), d);
        fmpz_sub(r, r, fmpz_mat_entry(cur, 2 * g - 1, g + k));
        fmpz_divexact(fmpz_mat_entry(u, g - 1, k), r, d);
        fmpz_set(fmpz_mat_entry(u, k, g - 1), fmpz_mat_entry(u, g - 1, k));
    }

    /* Todo: faster reduction using other entries of u as well? */
    sp2gz_trig(&vec[res], u);
    fmpz_mat_mul(cur, cur, &vec[res]);
    if (!fmpz_mat_is_one(&vec[res]))
    {
        res++;
    }

    /* Exchange c and d */
    sp2gz_j(&vec[res]);
    sp2gz_inv(&vec[res], &vec[res]);
    fmpz_mat_mul(next, cur, &vec[res]);
    res++;

    fmpz_mat_clear(cur);
    fmpz_mat_clear(u);
    fmpz_clear(d);
    fmpz_clear(r);
    return res;
}

static slong
sp2gz_get_parabolic(fmpz_mat_t next, fmpz_mat_struct* vec, const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    slong res = 0;
    fmpz_mat_t u, alpha;

    fmpz_mat_init(u, 2 * g, 2 * g);

    sp2gz_restrict(next, mat);
    sp2gz_embed(u, next);
    sp2gz_inv(u, u);
    fmpz_mat_mul(u, mat, u);

    fmpz_mat_window_init(alpha, u, 0, 0, g, g);
    if (!fmpz_mat_is_one(alpha))
    {
        sp2gz_block_diag(&vec[res], alpha);
        sp2gz_inv(&vec[res], &vec[res]);
        fmpz_mat_mul(u, u, &vec[res]);
        res++;
    }
    fmpz_mat_window_clear(alpha);

    fmpz_mat_window_init(alpha, u, 0, g, g, 2 * g);
    if (!fmpz_mat_is_zero(alpha))
    {
        sp2gz_trig(&vec[res], alpha);
        sp2gz_inv(&vec[res], &vec[res]);
        res++;
    }
    fmpz_mat_window_clear(alpha);

    fmpz_mat_clear(u);
    return res;
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
    }
    nb_max = 3 * nb_max + 3; /* for last delta reduction */

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

    /* Reduce delta one last time and recursive call */
    sp2gz_reduce_delta(cur, &gamma[nb_gamma], cur);
    if (!fmpz_mat_is_one(&gamma[nb_gamma]))
    {
        nb_gamma++;
    }
    if (g > 1)
    {
        fmpz_mat_init(w, 2 * (g - 1), 2 * (g - 1));
        add = sp2gz_get_parabolic(w, gamma + nb_gamma, cur);
        rec = sp2gz_decompose(&nb_rec, w);
        fmpz_mat_clear(w);
    }
    else
    {
        sp2gz_inv(&gamma[nb_gamma], cur);
        add = (fmpz_mat_is_one(&gamma[nb_gamma]) ? 0 : 1);
    }

    /* Make final vector */
    *nb = add + nb_rec + nb_gamma;
    res = flint_malloc(*nb * sizeof(fmpz_mat_struct));
    for (k = 0; k < *nb; k++)
    {
        fmpz_mat_init(&res[k], 2 * g, 2 * g);
    }

    for (k = 0; k < add; k++)
    {
        sp2gz_inv(&res[k], &gamma[nb_gamma + add - 1 - k]);
    }
    for (k = 0; k < nb_rec; k++)
    {
        sp2gz_embed(&res[add + k], &rec[k]);
    }
    for (k = 0; k < nb_gamma; k++)
    {
        sp2gz_inv(&res[add + nb_rec + k], &gamma[nb_gamma - 1 - k]);
    }

    fmpz_mat_clear(cur);
    for (k = 0; k < nb_max; k++)
    {
        fmpz_mat_clear(&gamma[k]);
    }
    flint_free(gamma);
    if (g > 1)
    {
        for (k = 0; k < nb_rec; k++)
        {
            fmpz_mat_clear(&rec[k]);
        }
        flint_free(rec);
    }
    return res;
}
