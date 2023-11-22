/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"
#include "acb_modular.h"
#include "acb_theta.h"

static slong
transform_kappa_g1(acb_t sqrtdet, const fmpz_mat_t mat, const fmpz_mat_t x,
    const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
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

    acb_mul_fmpz(sqrtdet, acb_mat_entry(tau, 0, 0), &y->c, prec);
    acb_add_fmpz(sqrtdet, sqrtdet, &y->d, prec);
    acb_sqrt(sqrtdet, sqrtdet, prec);

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
transform_kappa_j(acb_t sqrtdet, const fmpz_mat_t mat, const acb_mat_t tau, slong prec)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t gamma;
    acb_mat_t tau0;
    slong r, res;

    fmpz_mat_window_init(gamma, mat, g, 0, 2 * g, g);
    r = fmpz_mat_rank(gamma);
    fmpz_mat_window_clear(gamma);

    /* Mumford: theta_00(mtau) = det(tau0/i)^{1/2} theta_00(tau), and
       transform_sqrtdet(tau0) = i^{r/2} det(tau0/i)^{1/2} */
    acb_mat_window_init(tau0, tau, 0, 0, r, r);
    acb_theta_transform_sqrtdet(sqrtdet, tau0, prec);
    acb_mat_window_clear(tau0);

    res = -r;
    if (r % 2 == 1)
    {
        acb_mul_onei(sqrtdet, sqrtdet);
        res -= 2;
    }
    return res;
}

slong
acb_theta_transform_kappa(acb_t sqrtdet, const fmpz_mat_t mat,
    const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    fmpz_mat_struct * dec;
    fmpz_mat_t delta;
    fmpz_t det;
    slong nb_dec;
    fmpz_mat_t x;
    acb_mat_t w;
    acb_t c;
    slong k, res, e;
    ulong ab;

    fmpz_mat_init(x, 2, 2);
    acb_mat_init(w, g, g);
    acb_init(c);
    fmpz_init(det);
    dec = sp2gz_decompose(&nb_dec, mat);

    acb_one(sqrtdet);
    acb_mat_set(w, tau);
    res = 0;

    for (k = nb_dec - 1; k >= 0; k--)
    {
        if (sp2gz_is_trig(&dec[k]) || sp2gz_is_block_diag(&dec[k]))
        {
            /* theta_00(mtau) = theta_ab(tau) */
            fmpz_mat_window_init(delta, &dec[k], g, g, 2 * g, 2 * g);
            fmpz_mat_det(det, delta);
            fmpz_mat_window_clear(delta);

            if (fmpz_is_one(det))
            {
                acb_one(c);
            }
            else
            {
                acb_onei(c);
                res -= 2;
            }
        }
        else if (sp2gz_is_embedded(x, &dec[k]))
        {
            if (fmpz_cmp_si(fmpz_mat_entry(x, 1, 0), 0) < 0
                || (fmpz_is_zero(fmpz_mat_entry(x, 1, 0))
                    && fmpz_cmp_si(fmpz_mat_entry(x, 1, 1), 0) < 0))
            {
                fmpz_mat_neg(x, x);
                res += transform_kappa_g1(c, &dec[k], x, w, prec);
                acb_div_onei(c, c);
                res += 2;
            }
            else
            {
                res += transform_kappa_g1(c, &dec[k], x, w, prec);
            }
        }
        else /* embedded j */
        {
            res += transform_kappa_j(c, &dec[k], w, prec);
        }
        acb_siegel_transform(w, &dec[k], w, prec);
        acb_mul(sqrtdet, sqrtdet, c, prec);
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
    acb_mat_clear(w);
    acb_clear(c);
    for (k = 0; k < nb_dec; k++)
    {
        fmpz_mat_clear(&dec[k]);
    }
    flint_free(dec);
    return res & 7;
}
