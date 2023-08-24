/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

/* Use theta relation at low precision */

static void
acb_theta_transform_sqrtdet_lowprec(acb_t r, const acb_mat_t tau, const fmpz_mat_t mat)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong prec = ACB_THETA_LOW_PREC;
    slong kappa = acb_theta_transform_kappa(mat);
    slong b = -1;
    acb_mat_t w;
    acb_ptr z, th;
    acb_t t;
    fmpz_t eps;
    slong k;
    ulong ab;

    acb_mat_init(w, g, g);
    z = _acb_vec_init(g);
    th = _acb_vec_init(n);
    acb_init(t);
    fmpz_init(eps);

    while (b < 0)
    {
        /* Find b such that theta_{0,b}(gamma tau) is nonzero */
        prec *= 2;
        acb_mat_get_mid(w, tau);
        acb_siegel_transform(w, mat, w, prec);
        acb_theta_naive_0b(th, z, 1, w, prec);
        for (k = 0; k < n; k++)
        {
            if (!acb_contains_zero(&th[k]))
            {
                b = k;
                break;
            }
        }
    }

    ab = acb_theta_transform_char(eps, b, mat);
    acb_zero(r);

    while (acb_contains_zero(r))
    {
        acb_theta_naive_ind(th, b, z, 1, w, prec);
        acb_theta_naive_ind(r, ab, z, 1, tau, prec);
        acb_set_fmpz(t, eps);
        acb_add_si(t, t, kappa, prec);
        acb_mul_2exp_si(t, t, -2);
        acb_exp_pi_i(t, t, prec);
        acb_mul(r, r, t, prec);
        acb_div(r, &th[b], r, prec);
        prec *= 2;
    }

    acb_mat_clear(w);
    _acb_vec_clear(z, g);
    _acb_vec_clear(th, n);
    acb_clear(t);
    fmpz_clear(eps);
}

void
acb_theta_transform_sqrtdet(acb_t r, const fmpz_mat_t mat, const acb_mat_t tau, slong prec)
{
    acb_t x, y1, y2;

    acb_init(x);
    acb_init(y1);
    acb_init(y2);

    acb_siegel_cocycle_det(r, mat, tau, prec);
    acb_theta_transform_sqrtdet_lowprec(x, mat, tau);
    acb_sqrts(y1, y2, r, prec);

    if (acb_overlaps(y1, x) && acb_overlaps(y2, x))
    {
        arb_union(r, y1, y2, prec);
    }
    else if (acb_overlaps(y1, x))
    {
        arb_set(r, y1);
    }
    else if (acb_overlaps(y2, x))
    {
        arb_set(r, y2);
    }
    else
    {
        flint_printf("(acb_theta_transform_sqrtdet) Error: no overlap\n");
        flint_abort();
    }

    acb_clear(x);
    acb_clear(y1);
    acb_clear(y2);
}

/* static void */
/* acb_siegel_sqrtdet_propagate(acb_t r, acb_t r0, const acb_mat_t tau, const acb_mat_t tau0, */
/*     const fmpz_mat_t mat, slong prec, int lvl) */
/* { */
/*     slong g = acb_mat_nrows(tau); */
/*     acb_mat_t ctr, ball; */
/*     acb_t x; */
/*     slong j, k; */

/*     acb_mat_init(ctr, g, g); */
/*     acb_mat_init(ball, g, g); */
/*     acb_init(x); */

/*     flint_printf("Level %wd\n", lvl); */
/*     /\* flint_printf("Propagating between:\n"); */
/*     acb_mat_printd(tau0, 5); */
/*     acb_mat_printd(tau, 5); *\/ */

/*     /\* Set ctr to ball containing tau0, tau *\/ */
/*     acb_mat_add(ctr, tau, tau0, prec); */
/*     acb_mat_scalar_mul_2exp_si(ctr, ctr, -1); */
/*     acb_mat_set(ball, ctr); */
/*     for (j = 0; j < g; j++) */
/*     { */
/*         for (k = 0; k < g; k++) */
/*         { */
/*             acb_sub(x, acb_mat_entry(tau, j, k), acb_mat_entry(ball, j, k), prec); */
/*             arb_add_error(acb_realref(acb_mat_entry(ball, j, k)), acb_realref(x)); */
/*             arb_add_error(acb_imagref(acb_mat_entry(ball, j, k)), acb_imagref(x)); */
/*         } */
/*     } */
/*     if (!acb_mat_contains(ball, tau) || !acb_mat_contains(ball, tau0)) */
/*     { */
/*         flint_printf("(acb_siegel_cocycle_sqrtdet) Error: no matrix containment\n"); */
/*         acb_mat_printd(tau0, 5); */
/*         acb_mat_printd(tau, 5); */
/*         acb_mat_printd(ball, 5); */
/*         flint_abort(); */
/*     } */

/*     /\* Does det contain zero? If yes, recursion or fail, otherwise success *\/ */
/*     acb_siegel_cocycle_det(x, mat, ball, prec); */

/*     if (!acb_contains_zero(x)) */
/*     { */
/*         /\* Get the right square root, check that it contains r0 *\/ */
/*         acb_sqrt_no_cut(x, x, prec); */
/*         if (!acb_contains(x, r0)) */
/*         { */
/*             acb_neg(x, x); */
/*         } */
/*         if (!acb_contains(x, r0)) */
/*         { */
/*             flint_printf("(acb_siegel_cocycle_sqrtdet) Error: r0 not contained\n"); */
/*             acb_printd(x, 10); */
/*             flint_printf("\n"); */
/*             acb_printd(r0, 10); */
/*             flint_printf("\n"); */
/*             flint_abort(); */
/*         } */

/*         /\* Get the right square root det at tau *\/ */
/*         acb_siegel_cocycle_det(r, mat, tau, prec); */
/*         acb_sqrt_no_cut(r, r, prec); */
/*         if (!acb_contains(x, r)) */
/*         { */
/*             acb_neg(r, r); */
/*         } */
/*         if (!acb_contains(x, r)) */
/*         { */
/*             flint_printf("(acb_siegel_cocycle_sqrtdet) Error: r not contained\n"); */
/*             acb_printd(x, 10); */
/*             flint_printf("\n"); */
/*             acb_printd(r, 10); */
/*             flint_printf("\n"); */
/*             flint_abort(); */
/*         } */
/*     } */

/*     else if (acb_mat_overlaps(tau, tau0)) */
/*     { */
/*         acb_indeterminate(r); */
/*     } */

/*     else /\* x contains zero and no overlap: recursion *\/ */
/*     { */
/*         acb_siegel_sqrtdet_propagate(x, r0, ctr, tau0, mat, prec, lvl+1); */
/*         acb_siegel_sqrtdet_propagate(r, x, tau, ctr, mat, prec, lvl+1); */
/*     } */

/*     acb_mat_clear(ctr); */
/*     acb_mat_clear(ball); */
/*     acb_clear(x); */
/* } */

/* void */
/* acb_siegel_cocycle_sqrtdet(acb_t r, const fmpz_mat_t mat, const acb_mat_t tau, slong prec) */
/* { */
/*     slong g = acb_mat_nrows(tau); */
/*     acb_mat_t tau0; */
/*     acb_t x; */
/*     slong lp = ACB_THETA_LOW_PREC; */
/*     slong max = 16; */
/*     slong k; */

/*     acb_mat_init(tau0, g, g); */
/*     acb_init(x); */

/*     /\* Estimate necessary low precision; abort if too high *\/ */
/*     for (k = 0; k < max; k++) */
/*     { */
/*         acb_siegel_cocycle_det(x, mat, tau, lp); */
/*         lp *= 2; */
/*         if (!acb_contains_zero(x)) */
/*         { */
/*             flint_printf("Low prec = %wd\n", lp); */
/*             break; */
/*         } */
/*     } */

/*     if (k == max) */
/*     { */
/*         acb_indeterminate(r); */
/*     } */
/*     else */
/*     { */
/*         acb_mat_onei(tau0); */
/*         acb_siegel_sqrtdet_i(x, mat); */
/*         acb_siegel_sqrtdet_propagate(x, x, tau, tau0, mat, lp, 0); */

/*         acb_siegel_cocycle_det(r, mat, tau, prec); */
/*         acb_sqrt_no_cut(r, r, prec); */
/*         if (!acb_overlaps(x, r)) */
/*         { */
/*             acb_neg(r, r); */
/*         } */
/*         if (!acb_overlaps(x, r)) */
/*         { */
/*             flint_printf("(acb_siegel_cocycle_sqrtdet) Error: no overlap for final r\n"); */
/*             acb_printd(x, 10); */
/*             flint_printf("\n"); */
/*             acb_printd(r, 10); */
/*             flint_printf("\n"); */
/*             flint_abort(); */
/*         } */
/*     } */

/*     acb_mat_clear(tau0); */
/*     acb_clear(x); */
/* } */
