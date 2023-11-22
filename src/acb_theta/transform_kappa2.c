/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_modular.h"
#include "acb_theta.h"

static slong
transform_kappa2_g1(const fmpz_mat_t mat, const fmpz_mat_t x)
{
    slong g = sp2gz_dim(mat);
    psl2z_t y;
    int R[4];
    int S[4];
    int C;
    ulong ab;
    slong e, res;

    psl2z_init(y);

    /* set y to corresponding psl2z_t and use acb_modular_theta_transform */
    fmpz_set(&y->a, fmpz_mat_entry(x, 0, 0));
    fmpz_set(&y->b, fmpz_mat_entry(x, 0, 1));
    fmpz_set(&y->c, fmpz_mat_entry(x, 1, 0));
    fmpz_set(&y->d, fmpz_mat_entry(x, 1, 1));

    acb_modular_theta_transform(R, S, &C, y);

    /* find out where theta_00 is going */
    if (S[2] == 1) /* theta_2 */
    {
        ab = 1 << (2 * g - 1);
    }
    else if (S[2] == 2) /* theta_0 */
    {
        ab = 0;
    }
    else /* theta_1, since -theta_3 cannot happen (odd) */
    {
        ab = 1 << (g - 1);
    }
    acb_theta_transform_char(&e, mat, ab);

    /* adjust root of unity based on R */
    if (fmpz_is_zero(&y->c))
    {
        res = -R[2] - e;
    }
    else
    {
        res = -R[2] - 1 - e;
    }

    psl2z_clear(y);
    return res;
}

static slong
transform_kappa2_j(const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t gamma;
    slong r, res;

    fmpz_mat_window_init(gamma, mat, g, 0, 2 * g, g);
    r = fmpz_mat_rank(gamma);
    fmpz_mat_window_clear(gamma);

    res = -r;
    if (r % 2 == 1)
    {
        res -= 2;
    }
    return res;
}

slong
acb_theta_transform_kappa2(const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_struct * dec;
    fmpz_mat_t delta;
    fmpz_t det;
    slong nb_dec;
    fmpz_mat_t x;
    slong k, res, e;
    ulong ab;

    fmpz_mat_init(x, 2, 2);
    fmpz_init(det);
    dec = sp2gz_decompose(&nb_dec, mat);

    res = 0;
    for (k = nb_dec - 1; k >= 0; k--)
    {
        if (sp2gz_is_trig(&dec[k]) || sp2gz_is_block_diag(&dec[k]))
        {
            /* theta_00(mtau) = theta_ab(tau) */
            fmpz_mat_window_init(delta, &dec[k], g, g, 2 * g, 2 * g);
            fmpz_mat_det(det, delta);
            fmpz_mat_window_clear(delta);

            if (!fmpz_is_one(det))
            {
                res += 2;
            }
        }
        else if (sp2gz_is_embedded(x, &dec[k]))
        {
            if (fmpz_cmp_si(fmpz_mat_entry(x, 1, 0), 0) < 0
                || (fmpz_is_zero(fmpz_mat_entry(x, 1, 0))
                    && fmpz_cmp_si(fmpz_mat_entry(x, 1, 1), 0) < 0))
            {
                fmpz_mat_neg(x, x);
                res += transform_kappa2_g1(&dec[k], x);
                res += 2;
            }
            else
            {
                res += transform_kappa2_g1(&dec[k], x);
            }
        }
        else /* embedded j */
        {
            res += transform_kappa2_j(&dec[k]);
        }
    }

    /* Adjust final sign based on transformation of coordinates */
    acb_theta_transform_char(&e, mat, 0);
    res -= e;
    ab = 0;
    for (k = 0; k < nb_dec; k++)
    {
        ab = acb_theta_transform_char(&e, &dec[k], ab);
        res += e;
    }

    fmpz_mat_clear(x);
    for (k = 0; k < nb_dec; k++)
    {
        fmpz_mat_clear(&dec[k]);
    }
    flint_free(dec);
    return res & 3;
}

