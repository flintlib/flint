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
fmpz_mat_bound_inf_norm(mag_t b, const fmpz_mat_t mat)
{
    slong r = acb_mat_nrows(mat);
    slong c = acb_mat_ncols(mat);
    arb_mat_t m;

    arb_mat_init(m, r, c);
    arb_mat_set_fmpz_mat(m, mat);
    arb_mat_bound_inf_norm(b, m);
    arb_mat_clear(m);
}

/* Todo: g * (...) is an emergency fix, what is the right value here? */
static slong
acb_siegel_reduce_real_lowprec(const mag_t ntau, const mag_t nmat, slong g, slong prec)
{
    slong lp = ACB_THETA_LOW_PREC;
    slong res;
    mag_t b;

    mag_init(b);
    mag_mul(b, ntau, nmat);
    res = FLINT_MIN(prec, g * (lp + FLINT_MAX(0, mag_get_d_log2_approx(b))));
    mag_clear(b);

    return res;
}

static slong
acb_siegel_reduce_imag_lowprec(const mag_t ntau, const mag_t ndet, const mag_t nmat,
    slong g, slong prec)
{
    slong lp = ACB_THETA_LOW_PREC;
    slong res;
    mag_t b;

    mag_init(b);
    mag_mul(b, ntau, nmat);
    mag_mul(b, b, b);
    mag_mul(b, b, ntau);
    mag_div(b, b, ndet);
    res = FLINT_MIN(prec, g * (lp + FLINT_MAX(0, mag_get_d_log2_approx(b))));
    mag_clear(b);

    return res;
}

void
acb_siegel_reduce(acb_mat_t res, fmpz_mat_t mat, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong lp;
    fmpz_mat_t m;
    acb_mat_t cur;
    arb_mat_t im;
    acb_t det;
    arb_t abs;
    arb_t t;
    mag_t ntau, nmat, ndet;

    int stop = 0;
    slong j, j0;

    fmpz_mat_init(m, 2 * g, 2 * g);
    acb_mat_init(cur, g, g);
    arb_mat_init(im, g, g);
    acb_init(det);
    arb_init(abs);
    arb_init(t);
    mag_init(ntau);
    mag_init(nmat);
    mag_init(ndet);

    acb_mat_bound_inf_norm(ntau, tau);
    acb_mat_get_imag(im, tau);
    arb_mat_det(abs, im, prec);
    arb_get_mag_lower(ndet, abs);
    if (mag_is_inf(ntau) || mag_is_zero(ndet))
    {
        stop = 1;
    }

    fmpz_mat_one(mat);
    while (!stop)
    {
        /* Choose precision, reduce imaginary part */
        fmpz_mat_bound_inf_norm(nmat, mat);
        lp = acb_siegel_reduce_imag_lowprec(ntau, ndet, nmat, g, prec);
        acb_siegel_transform(cur, mat, tau, lp);
        acb_siegel_reduce_imag(m, cur, lp);
        fmpz_mat_mul(mat, m, mat);

        /* Choose precision, reduce real part */
        fmpz_mat_bound_inf_norm(nmat, mat);
        lp = acb_siegel_reduce_real_lowprec(ntau, nmat, g, prec);
        acb_siegel_transform(cur, m, cur, lp);
        acb_siegel_reduce_real(m, cur, lp);
        fmpz_mat_mul(mat, m, mat);

        /* Loop over fundamental matrices (keeping same precision) */
        acb_siegel_transform(cur, m, cur, lp);
        j0 = -1;
        arb_one(t);
        for (j = 0; j < sp2gz_nb_fundamental(g); j++)
        {
            sp2gz_fundamental(m, j);
            acb_siegel_cocycle_det(det, m, cur, lp);
            acb_abs(abs, det, lp);
            if (arb_lt(abs, t))
            {
                j0 = j;
                arb_set(t, abs);
            }
        }

        /* Apply fundamental matrix if found */
        if (j0 != -1)
        {
            sp2gz_fundamental(m, j0);
            fmpz_mat_mul(mat, m, mat);
        }
        else
        {
            stop = 1;
        }
    }

    /* Final transform at full precision */
    acb_siegel_transform(res, mat, tau, prec);

    fmpz_mat_clear(m);
    acb_mat_clear(cur);
    arb_mat_clear(im);
    acb_clear(det);
    arb_clear(abs);
    arb_clear(t);
    mag_clear(ntau);
    mag_clear(nmat);
    mag_clear(ndet);
}
