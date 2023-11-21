/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

#define ACB_SIEGEL_REDUCE_MAG_BOUND 1000000

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

/* Todo: better choice of precision here? */
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

static void
acb_siegel_reduce_real(fmpz_mat_t mat, const acb_mat_t tau)
{
    slong g = acb_mat_nrows(tau);
    slong j, k;
    fmpz_t c;

    fmpz_init(c);

    fmpz_mat_one(mat);
    for (j = 0; j < g; j++)
    {
        for (k = j; k < g; k++)
        {
            /* this must succeed given the bounds on ndet and ntau */
            arf_get_fmpz(c, arb_midref(acb_realref(acb_mat_entry(tau, j, k))),
                         ARF_RND_NEAR);
            fmpz_neg(fmpz_mat_entry(mat, j, k + g), c);
        }
        for (k = 0; k < j; k++)
        {
            fmpz_set(fmpz_mat_entry(mat, j, k + g),
                     fmpz_mat_entry(mat, k, j + g));
        }
    }

    fmpz_clear(c);
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

static void
acb_siegel_reduce_imag(fmpz_mat_t mat, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_mat_t im;
    fmpz_mat_t U;

    arb_mat_init(im, g, g);
    fmpz_mat_init(U, g, g);

    acb_mat_get_imag(im, tau);
    arb_mat_spd_lll_reduce(U, im, prec);
    sp2gz_block_diag(mat, U);

    arb_mat_clear(im);
    fmpz_mat_clear(U);
}

void
acb_siegel_reduce(fmpz_mat_t mat, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong lp;
    fmpz_mat_t m;
    acb_mat_t w, c;
    arb_mat_t im;
    acb_t det;
    arb_t abs;
    arb_t t;
    mag_t ntau, nmat, ndet;

    int stop = 0;
    slong j, j0;

    fmpz_mat_init(m, 2 * g, 2 * g);
    acb_mat_init(w, g, g);
    acb_mat_init(c, g, g);
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
    if (mag_cmp_2exp_si(ntau, ACB_SIEGEL_REDUCE_MAG_BOUND) >= 0
        || mag_cmp_2exp_si(ndet, -ACB_SIEGEL_REDUCE_MAG_BOUND) <= 0)
    {
        stop = 1;
    }

    fmpz_mat_one(mat);
    while (!stop)
    {
        /* Choose precision, reduce imaginary part */
        fmpz_mat_bound_inf_norm(nmat, mat);
        lp = acb_siegel_reduce_imag_lowprec(ntau, ndet, nmat, g, prec);
        acb_siegel_transform(w, mat, tau, lp);
        acb_siegel_reduce_imag(m, w, lp);
        fmpz_mat_mul(mat, m, mat);

        /* Choose precision, check transform is reduced, reduce real part */
        fmpz_mat_bound_inf_norm(nmat, mat);
        lp = acb_siegel_reduce_real_lowprec(ntau, nmat, g, prec);
        acb_siegel_transform(w, m, w, lp);
        acb_mat_get_imag(im, w);
        if (!arb_mat_spd_is_lll_reduced(im, -10, lp))
        {
            stop = 1;
            break;
        }
        acb_siegel_reduce_real(m, w);
        fmpz_mat_mul(mat, m, mat);

        /* Loop over fundamental matrices (keeping same precision) */
        acb_siegel_transform(w, m, w, lp);
        j0 = -1;
        arb_one(t);
        for (j = 0; j < sp2gz_nb_fundamental(g); j++)
        {
            sp2gz_fundamental(m, j);
            acb_siegel_cocycle(c, m, w, lp);
            acb_mat_det(det, c, lp);
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

    fmpz_mat_clear(m);
    acb_mat_clear(w);
    acb_mat_clear(c);
    arb_mat_clear(im);
    acb_clear(det);
    arb_clear(abs);
    arb_clear(t);
    mag_clear(ntau);
    mag_clear(nmat);
    mag_clear(ndet);
}
