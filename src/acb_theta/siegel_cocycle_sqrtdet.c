/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

/* todo: move out? also used in agm_sqrt */
static void
acb_sqrt_no_cut(acb_t r, const acb_t x, slong prec)
{
    if (arb_contains_zero(acb_imagref(x)) && arb_is_negative(acb_realref(x)))
    {
        acb_neg(r, x);
        acb_sqrt(r, r, prec);
        acb_mul_onei(r, r);
    }
    else
    {
        acb_sqrt(r, x, prec);
    }
}

static void
acb_siegel_sqrtdet_i(acb_t r, const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    slong prec = ACB_THETA_LOW_PREC;
    acb_mat_t tau;

    acb_mat_init(tau, g, g);

    acb_mat_onei(tau);
    acb_siegel_cocycle_det(r, mat, tau, prec);

    /* r should be exact */
    if (!acb_is_exact(r))
    {
        flint_printf("(acb_siegel_cocycle_sqrtdet) Error: not exact\n");
        acb_printd(r, 5);
        flint_printf("\n");
        flint_abort();
    }
    acb_sqrt_no_cut(r, r, prec);

    acb_mat_clear(tau);
}

static void
acb_siegel_sqrtdet_propagate(acb_t r, acb_t r0, const acb_mat_t tau, const acb_mat_t tau0,
    const fmpz_mat_t mat, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_mat_t ctr, ball;
    acb_t x;
    slong j, k;

    acb_mat_init(ctr, g, g);
    acb_mat_init(ball, g, g);
    acb_init(x);

    /* Set ctr to ball containing tau0, tau */
    acb_mat_add(ctr, tau, tau0, prec);
    acb_mat_scalar_mul_2exp_si(ctr, ctr, -1);
    acb_mat_set(ball, ctr);
    for (j = 0; j < g; j++)
    {
        for (k = 0; k < g; k++)
        {
            acb_sub(x, acb_mat_entry(tau, j, k), acb_mat_entry(ball, j, k), prec);
            arb_add_error(acb_realref(acb_mat_entry(ball, j, k)), acb_realref(x));
            arb_add_error(acb_imagref(acb_mat_entry(ball, j, k)), acb_imagref(x));
        }
    }
    if (!acb_mat_contains(ball, tau) || !acb_mat_contains(ball, tau0))
    {
        flint_printf("(acb_siegel_cocycle_sqrtdet) Error: no matrix containment\n");
        acb_mat_printd(tau0, 5);
        acb_mat_printd(tau, 5);
        acb_mat_printd(ball, 5);
        flint_abort();
    }

    /* Does det contain zero? If yes, recursion or fail, otherwise success */
    acb_siegel_cocycle_det(x, mat, ball, prec);

    if (!acb_contains_zero(x))
    {
        /* Get the right square root, check that it contains r0 */
        acb_sqrt_no_cut(x, x, prec);
        if (!acb_contains(x, r0))
        {
            acb_neg(x, x);
        }
        if (!acb_contains(x, r0))
        {
            flint_printf("(acb_siegel_cocycle_sqrtdet) Error: r0 not contained\n");
            acb_printd(x, 10);
            flint_printf("\n");
            acb_printd(r0, 10);
            flint_printf("\n");
            flint_abort();
        }

        /* Get the right square root det at tau */
        acb_siegel_cocycle_det(r, mat, tau, prec);
        acb_sqrt_no_cut(r, r, prec);
        if (!acb_contains(x, r))
        {
            acb_neg(r, r);
        }
        if (!acb_contains(x, r))
        {
            flint_printf("(acb_siegel_cocycle_sqrtdet) Error: r not contained\n");
            acb_printd(x, 10);
            flint_printf("\n");
            acb_printd(r, 10);
            flint_printf("\n");
            flint_abort();
        }
    }

    else if (acb_mat_overlaps(tau, tau0))
    {
        acb_indeterminate(r);
    }

    else /* x contains zero and no overlap: recursion */
    {
        acb_siegel_sqrtdet_propagate(x, r0, ctr, tau0, mat, prec);
        acb_siegel_sqrtdet_propagate(r, x, tau, ctr, mat, prec);
    }

    acb_mat_clear(ctr);
    acb_mat_clear(ball);
    acb_clear(x);
}

void
acb_siegel_cocycle_sqrtdet(acb_t r, const fmpz_mat_t mat, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_mat_t tau0;
    acb_t x;
    slong lp = ACB_THETA_LOW_PREC;
    slong max = 16;
    slong k;

    acb_mat_init(tau0, g, g);
    acb_init(x);

    /* Estimate necessary low precision; abort if too high */
    for (k = 0; k < max; k++)
    {
        acb_siegel_cocycle_det(x, mat, tau, lp);
        lp *= 2;
        if (!acb_contains_zero(x))
        {
            break;
        }
    }

    if (k == max)
    {
        acb_indeterminate(r);
    }
    else
    {
        acb_mat_onei(tau0);
        acb_siegel_sqrtdet_i(x, mat);
        acb_siegel_sqrtdet_propagate(x, x, tau, tau0, mat, lp);

        acb_siegel_cocycle_det(r, mat, tau, prec);
        acb_sqrt_no_cut(r, r, prec);
        if (!acb_contains(x, r))
        {
            acb_neg(r, r);
        }
        if (!acb_contains(x, r))
        {
            flint_printf("(acb_siegel_cocycle_sqrtdet) Error: final r not contained\n");
            acb_printd(x, 10);
            flint_printf("\n");
            acb_printd(r, 10);
            flint_printf("\n");
            flint_abort();
        }
    }

    acb_mat_clear(tau0);
    acb_clear(x);
}
