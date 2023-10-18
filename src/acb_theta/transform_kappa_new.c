/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static slong
transform_kappa_g1(acb_t sqrtdet, const fmpz_mat_t mat, const fmpz_mat_t x,
    const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_mat_t c;
    psl2z_t y;
    int R[4];
    int S[4];
    int C;
    ulong ab;
    slong e, res;

    acb_mat_init(c, g, g);
    psl2z_init(y);

    /* set y to corresponding psl2z_t and use acb_modular_theta_transform */
    if (fmpz_cmp_si(fmpz_mat_entry(x, 1, 0), 0) < 0
        || (fmpz_is_zero(fmpz_mat_entry(x, 1, 0))
            && fmpz_cmp_si(fmpz_mat_entry(x, 1, 1), 0) < 0))
    {
        fmpz_mat_neg(x, x);
    }
    fmpz_set(&y->a, fmpz_mat_entry(x, 0, 0));
    fmpz_set(&y->b, fmpz_mat_entry(x, 0, 1));
    fmpz_set(&y->c, fmpz_mat_entry(x, 1, 0));
    fmpz_set(&y->d, fmpz_mat_entry(x, 1, 1));

    acb_modular_theta_transform(&R, &S, &C, y);

    acb_siegel_cocycle(c, mat, tau, prec);
    acb_mat_det(sqrtdet, c, prec);
    acb_div_onei(sqrtdet, sqrtdet); /* lies in right half plane */
    acb_sqrt(det, det);
    acb_mul(sqrtdet, det, prec);

    /* find out where theta_00 is going */
    if (S[2] == 0) /* -theta_3 */
    {
        ab = (1 << (2 * g - 1)) + (1 << (g - 1));
        res += 4;
    }
    else if (S[2] == 1) /* theta_2 */
    {
        ab = 1 << (2 * g - 1);
    }
    else if (S[2] == 2) /* theta_0 */
    {
        ab = 0;
    }
    else /* theta_1 */
    {
        ab = 1 << (g - 1);
    }
    ab = acb_theta_transform_char(&e, mat, ab);
    if (ab != 0)
    {
        flint_printf("(transform_kappa_new) error g1: expected ab = 0, got %wd\n", ab);
    }

    /* adjust root of 1 based on R */
    res = -R[2] - e;

    acb_mat_clear(c);
    psl2z_clear(y);
}

static slong
transform_kappa_j(acb_t sqrtdet, const fmpz_mat_t mat, const acb_mat_t tau, slong prec)
{
    
}


slong acb_theta_transform_kappa_new(acb_t sqrtdet, const fmpz_mat_t mat,
    const acb_mat_t tau, slong prec)
{
    slong nb_dec;
    fmpz_mat_struct* dec;
    fmpz_mat_t x;
    acb_mat_t w, c, tau0;
    acb_t det;
    slong k, e, res;
    ulong ab;

    fmpz_mat_init(x, 2, 2);
    acb_mat_init(w, g, g);
    acb_mat_init(c, g, g);
    acb_init(det);
    dec = sp2gz_decompose(&nb_dec, mat);

    acb_one(sqrtdet);
    acb_set(w, tau);

    for (k = nb_dec - 1; k >= 0; k--)
    {
        if (sp2gz_is_trig(&dec[k]) || sp2gz_is_block_diag(&dec[k]))
        {
            /* theta_00(mtau) = theta_00(tau) */
            ab = acb_theta_transform_char(&e, &dec[k], 0);
            if (ab != 0)
            {
                flint_printf("(transform_kappa_new) error trig: expected ab = 0, got %wd\n", ab);
            }
            res -= e;
        }
        else if (sp2gz_is_embedded(x, &dec[k]))
        {
            res += transform_kappa_g1(det, &dec[k], x, w, prec);
            acb_mul(sqrtdet, sqrtdet, det, prec);
        }
        else /* embedded j */
        {
            res += transform_kappa_j(det, &dec[k], w, prec);
            acb_mul(sqrtdet, sqrtdet, det, prec);
        }
        acb_siegel_transform(w, &dec[k], w, prec);
    }
}
