/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static int
accept_naive(const acb_mat_t tau, slong prec)
{
    arb_t test;
    int res;

    arb_init(test);

    arb_set(test, acb_imagref(acb_mat_entry(tau, 0, 0)));
    arb_mul_si(test, test, ACB_THETA_NAIVE_THRESHOLD, prec);
    arb_sub_si(test, test, prec, prec);

    res = arb_is_positive(test);

    arb_clear(test);
    return res;
}

static int
is_reduced_z(acb_srcptr z, slong g, slong prec)
{
    arb_t abs, m;
    slong k;
    int res;

    arb_init(abs);
    arb_init(m);

    arb_zero(m);
    for (k = 0; k < g; k++)
    {
        acb_abs(abs, &z[k], prec);
        arb_max(m, m, abs, prec);
    }

    arb_mul_si(m, m, ACB_THETA_REDUCE_Z, prec);
    arb_sub_si(m, m, 1, prec);
    res = !arb_is_positive(m);

    arb_clear(abs);
    arb_clear(m);
    return res;
}


static void
theta2_naive(acb_ptr th2, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_ptr z_ext;

    z_ext = _acb_vec_init(2 * g);
    _acb_vec_set(z_ext, z, g);

    acb_theta_naive_all(th2, z_ext, 2, tau, prec);
    acb_theta_vecsqr(th2, th2, 1 << (2 * g + 1), prec);

    _acb_vec_clear(z_ext, 2 * g);
}

static void
theta2_newton_dupl_z(acb_ptr th2, acb_srcptr z, const acb_mat_t tau,
                     slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << (2 * g + 1);
    acb_ptr zmod;
    acb_ptr roots;
    slong k;
    slong cnt = 0;
    slong lowprec = acb_theta_balance_lowprec(g, prec);
    int r;

    zmod = _acb_vec_init(2 * g);
    roots = _acb_vec_init(n);

    /* Reduce z */
    _acb_vec_set(zmod, z, g);
    while (!is_reduced_z(zmod, g, prec))
    {
        _acb_vec_scalar_mul_2exp_si(zmod, zmod, g, -1);
        cnt++;
    }

    /* Attempt to collect nonzero, lowprec square roots */
    acb_theta_naive_all(roots, zmod, 2, tau, lowprec);

    flint_printf("Reduced z:\n");
    for (k = 0; k < g; k++)
    {
        acb_printd(&zmod[k], 10);
        flint_printf("\n");
    }
    flint_printf("Lowprec theta values for reduced z:\n");
    for (k = 0; k < n; k++)
    {
        acb_printd(&roots[k], 10);
        flint_printf("\n");
    }

    /* If one of the relevant square roots are zero, fall back to naive */
    r = 1;
    for (k = 1 << (2 * g); k < (1 << (2 * g)) + (1 << g); k++)
    {
        if (acb_contains_zero(&roots[k]))
            r = 0;
    }
    for (k = 1 << (2 * g); k < 1 << (2 * g + 1); k += (1 << g))
    {
        if (acb_contains_zero(&roots[k]))
            r = 0;
    }
    for (k = 0; k < 1 << (2 * g); k++)
    {
        if (acb_contains_zero(&roots[k]))
            r = 0;
    }

    if (!r)
    {
        theta2_naive(th2, z, tau, prec);
    }
    else
    {
        acb_theta_newton_all_sqr(th2, zmod, tau, prec);

        flint_printf("Highprec values for reduced z:\n");
        for (k = 0; k < n; k++)
        {
            acb_printd(&th2[k], 10);
            flint_printf("\n");
        }
        for (k = 0; k < n; k++)
        {
            acb_theta_agm_sqrt_lowprec(&th2[k], &th2[k], &roots[k], prec);
        }
        for (k = 0; k < cnt; k++)
        {
            acb_theta_dupl_z(th2, th2, g, prec);
        }
        acb_theta_vecsqr(th2, th2, n, prec);
    }

    _acb_vec_clear(zmod, 2 * g);
    _acb_vec_clear(roots, n);
}

static void
theta2_newton(acb_ptr th2, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_ptr new_z;
    arb_mat_t v;
    arb_mat_t Y;
    acb_mat_t col, row, prod;
    acb_t c;
    slong *r1;
    slong *r2;
    slong k;

    new_z = _acb_vec_init(g);
    arb_mat_init(v, g, 1);
    arb_mat_init(Y, g, g);
    acb_mat_init(col, g, 1);
    acb_mat_init(row, 1, g);
    acb_mat_init(prod, 1, 1);
    acb_init(c);
    r1 = flint_malloc(g * sizeof(slong));
    r2 = flint_malloc(g * sizeof(slong));

    acb_mat_get_imag(Y, tau);
    arb_mat_inv(Y, Y, prec);

    /* Approximate imaginary part */
    for (k = 0; k < g; k++)
    {
        arb_set(arb_mat_entry(v, k, 0), acb_imagref(&z[k]));
    }
    arb_mat_mul(v, Y, v, prec);
    acb_theta_eld_round(r1, v);

    /* Substract from z */
    for (k = 0; k < g; k++)
    {
        acb_set_si(arb_mat_entry(col, k, 0), r1[k]);
    }
    acb_mat_mul(col, tau, col, prec);
    for (k = 0; k < g; k++)
    {
        acb_sub(&new_z[k], &z[k], acb_mat_entry(col, k, 0), prec);
    }

    /* Approximate real part and substract */
    for (k = 0; k < g; k++)
    {
        arb_set(arb_mat_entry(v, k, 0), acb_realref(&new_z[k]));
    }
    acb_theta_eld_round(r2, v);
    for (k = 0; k < g; k++)
    {
        acb_sub_si(&new_z[k], &new_z[k], r2[k], prec);
    }

    flint_printf("(all_sqr) Newton: z reduction is\n");
    for (k = 0; k < g; k++)
    {
        acb_printd(&new_z[k], 10);
        flint_printf("\n");
    }
    flint_printf("r1: %wd", r1[0]);
    for (k = 1; k < g; k++)
    {
        flint_printf(", %wd", r1[k]);
    }
    flint_printf("\nr2: %wd", r2[0]);
    for (k = 1; k < g; k++)
    {
        flint_printf(", %wd", r2[k]);
    }
    flint_printf("\n");

    /* Get theta from dupl_z */
    theta2_newton_dupl_z(th2, new_z, tau, prec);

    flint_printf("(all_sqr) Newton: before exponentials\n");
    for (k = 0; k < (1 << (2 * g + 1)); k++)
    {
        acb_printd(&th2[k], 10);
        flint_printf("\n");
    }

    /* Get multiplicative factor */
    for (k = 0; k < g; k++)
    {
        acb_set_si(acb_mat_entry(col, k, 0), r1[k]);
    }
    acb_mat_transpose(row, col);
    acb_mat_mul(col, tau, col, prec);
    acb_mat_mul(prod, row, col, prec);
    acb_set(c, acb_mat_entry(prod, 0, 0));

    for (k = 0; k < g; k++)
    {
        acb_set(acb_mat_entry(col, k, 0), &new_z[k]);
        acb_set_si(acb_mat_entry(row, 0, k), r1[k]);
    }
    acb_mat_mul(prod, row, col, prec);
    acb_addmul_si(c, acb_mat_entry(prod, 0, 0), 2, prec);

    acb_neg(c, c);
    acb_mul_2exp_si(c, c, 1);
    acb_exp_pi_i(c, c, prec);
    _acb_vec_scalar_mul(th2, th2, 1 << (2 * g), c, prec);

    flint_printf("(all_sqr) Newton: after exponentials\n");
    for (k = 0; k < (1 << (2 * g + 1)); k++)
    {
        acb_printd(&th2[k], 10);
        flint_printf("\n");
    }

    _acb_vec_clear(new_z, g);
    arb_mat_clear(v);
    arb_mat_clear(Y);
    acb_mat_clear(col);
    acb_mat_clear(row);
    acb_mat_clear(prod);
    acb_clear(c);
    flint_free(r1);
    flint_free(r2);
}

static int
lowprec_roots(acb_ptr roots, acb_srcptr z, const acb_mat_t tau,
              const fmpz_mat_t mat, slong lowprec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_ptr th;
    ulong k;
    int res = 1;

    th = _acb_vec_init(2 * n * n);

    acb_theta_naive_all(th, z, 2, tau, lowprec);
    acb_theta_transform_proj(roots, th, mat, lowprec);
    acb_theta_transform_proj(roots + n, th + n * n, mat, lowprec);

    flint_printf("Lowprec roots:\n");
    for (k = 0; k < 2 * n; k++)
    {
        acb_printd(&roots[k], 10);
        flint_printf("\n");
    }

    for (k = 0; k < n; k++)
    {
        if (acb_contains_zero(&roots[k]))
        {
            res = 0;
            break;
        }
    }

    _acb_vec_clear(th, 2 * n * n);
    return res;
}

void
acb_theta_all_sqr(acb_ptr th2, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    fmpz_mat_t D, R;
    acb_mat_t aux;
    acb_ptr zmod;
    acb_ptr roots;
    fmpz_t den;
    int res;
    slong j0, k;
    slong lowprec = acb_theta_balance_lowprec(g, prec);

    fmpz_mat_init(D, 2 * g, 2 * g);
    fmpz_mat_init(R, 2 * g, 2 * g);
    acb_mat_init(aux, g, g);
    zmod = _acb_vec_init(2 * g);
    roots = _acb_vec_init(2 * n);
    fmpz_init(den);

    res = accept_naive(tau, prec);
    if (res)
    {
        flint_printf("(all_sqr) Fall back to naive.\n");
        theta2_naive(th2, z, tau, prec);
        goto exit;
    }

    res = acb_theta_is_balanced(&j0, tau, prec);
    if (res)
    {
        flint_printf("(all_sqr) Fall back to newton.\n");
        theta2_newton(th2, z, tau, prec);
        goto exit;
    }

    acb_theta_balance(zmod, aux, D, z, tau, j0);

    flint_printf("Before real reduction:\n");
    acb_mat_printd(aux, 10);
    for (k = 0; k < g; k++)
    {
        acb_printd(&zmod[k], 10);
        flint_printf("\n");
    }

    /* Reduce real part */
    acb_siegel_reduce_real(R, aux, prec);
    acb_siegel_transform(aux, R, aux, prec);

    flint_printf("After real reduction:\n");
    acb_mat_printd(aux, 10);
    for (k = 0; k < g; k++)
    {
        acb_printd(&zmod[k], 10);
        flint_printf("\n");
    }

    /* Compute th2 recursively, act by R^-1 */
    acb_theta_all_sqr(th2, zmod, aux, prec);

    flint_printf("At exit of recursive call:\n");
    for (k = 0; k < 2 * n * n; k++)
    {
        acb_printd(&th2[k], 10);
        flint_printf("\n");
    }
    flint_printf("\n");

    fmpz_mat_inv(R, den, R);
    acb_theta_transform_all_sqr_proj(th2, th2, R, prec);
    acb_theta_transform_all_sqr_proj(th2 + n * n, th2 + n * n, R, prec);
    acb_siegel_transform(aux, R, aux, prec);

    flint_printf("After reverse real reduction:\n");
    for (k = 0; k < 2 * n * n; k++)
    {
        acb_printd(&th2[k], 10);
        flint_printf("\n");
    }
    flint_printf("\n");

    /* Attempt duplication at D.aux: get roots at low precision */
    res = lowprec_roots(roots, zmod, aux, D, lowprec);

    if (!res)
    {
        /* Some roots are zero: abort with naive algorithm */
        theta2_naive(th2, z, tau, prec);
        goto exit;
    }

    /* Duplicate */
    acb_theta_transform_sqr_proj(th2, th2, D, prec);
    acb_theta_transform_sqr_proj(th2 + n, th2 + n * n, D, prec);

    flint_printf("Before duplication:\n");
    for (k = 0; k < 2 * n; k++)
    {
        acb_printd(&th2[k], 10);
        flint_printf("\n");
    }
    flint_printf("\n");

    flint_printf("Roots:\n");
    for (k = 0; k < 2 * n; k++)
    {
        acb_printd(&roots[k], 10);
        flint_printf("\n");
    }
    flint_printf("\n");

    for (k = 0; k < 2 * n; k++)
    {
        acb_theta_agm_sqrt_lowprec(&th2[k], &th2[k], &roots[k], prec);
    }
    acb_theta_dupl_all(th2, th2, g, prec);

    flint_printf("After duplication:\n");
    for (k = 0; k < n * n; k++)
    {
        acb_printd(&th2[k], 10);
        flint_printf("\n");
    }
    flint_printf("\n");

    /* Act by inverse of D and x2 */
    fmpz_mat_inv(D, den, D);
    acb_theta_transform_all_sqr_proj(th2, th2, D, prec);
    acb_theta_transform_all_sqr_proj(th2 + n * n, th2 + n * n, D, prec);
    _acb_vec_scalar_mul_2exp_si(th2, th2, 2 * n * n, j0 + 1);

    flint_printf("After final transformation:\n");
    for (k = 0; k < 2 * n * n; k++)
    {
        acb_printd(&th2[k], 10);
        flint_printf("\n");
    }
    flint_printf("\n");

    goto exit;

  exit:
    {
        fmpz_mat_clear(D);
        fmpz_mat_clear(R);
        acb_mat_clear(aux);
        _acb_vec_clear(zmod, g);
        _acb_vec_clear(roots, 2 * n);
        fmpz_clear(den);
    }
}
