/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

/* todo: move out */
void
acb_mat_inf_norm(arb_t res, const acb_mat_t A, slong prec)
{
    slong i, j, m, n;
    arb_t s, t;

    m = acb_mat_nrows(A);
    n = acb_mat_nrows(A);

    if (m == 0 || n == 0)
    {
        arb_zero(res);
        return;
    }

    arb_init(s);
    arb_init(t);

    arb_zero(res);

    for (i = 0; i < m; i++)
    {
        acb_abs(s, acb_mat_entry(A, i, 0), prec);

        for (j = 1; j < n; j++)
        {
            acb_abs(t, acb_mat_entry(A, i, j), prec);
            arb_add(s, s, t, prec);
        }

        arb_max(res, res, s, prec);
    }

    arb_clear(s);
    arb_clear(t);
}

static void
diagonal_certify(arb_t epsilon, arb_t eta1, arb_t eta2, const acb_mat_t D, const acb_mat_t H, slong prec)
{
    arb_t mu, sigma, alpha, t, u, v;
    acb_t d;
    slong i, j, n;

    arb_init(mu);
    arb_init(sigma);
    arb_init(alpha);
    arb_init(t);
    arb_init(u);
    arb_init(v);
    acb_init(d);

    n = acb_mat_nrows(D);

    /* D = diagonal matrix; H = off diagonal matrix */

    /* mu = ||D|| */
    acb_mat_inf_norm(mu, D, prec);

    /* sigma = sigma(D) = separation number */
    arb_pos_inf(sigma);
    for (i = 0; i < n; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            acb_sub(d, acb_mat_entry(D, i, i), acb_mat_entry(D, j, j), prec);
            acb_abs(t, d, prec);
            arb_min(sigma, sigma, t, prec);
        }
    }

    /* eta1 = ||Delta(H)||   Delta = diagonal projection  */
    arb_zero(eta1);

    /* eta2 = ||Omega(H)||   Omega = off-diagonal projection  */
    acb_mat_inf_norm(eta2, H, prec);

    /* alpha = min(sigma / (6 mu), 1/4) */
    arb_div(t, sigma, mu, prec);
    arb_div_ui(t, t, 6, prec);
    arb_set_d(u, 0.25);
    arb_min(alpha, t, u, prec);

    arb_add(t, eta1, eta2, prec);
    arb_mul(u, alpha, mu, prec);
    arb_mul_2exp_si(u, u, -3);

    arb_mul(v, alpha, sigma, prec);
    arb_mul_2exp_si(v, v, -3);

    if (arb_le(t, u) && arb_le(eta2, v))
    {
        arb_div(epsilon, eta2, sigma, prec);
        arb_mul_ui(epsilon, epsilon, 3, prec);
    }
    else
    {
        arb_indeterminate(epsilon);
    }

    arb_clear(mu);
    arb_clear(sigma);
    arb_clear(alpha);
    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
    acb_clear(d);
}

int
acb_mat_eig_simple_vdhoeven_mourrain(acb_ptr E,
    acb_mat_t L, acb_mat_t R, const acb_mat_t A,
    acb_srcptr E_approx, const acb_mat_t R_approx, slong prec)
{
    acb_mat_t D, T, AT;
    int result;
    slong i, j, n;

    result = 0;
    n = acb_mat_nrows(A);

    if (n == 0)
        return 1;

    if (n == 1)
    {
        acb_set_round(E, acb_mat_entry(A, 0, 0), prec);
        if (L != NULL)
            acb_one(acb_mat_entry(L, 0, 0));
        if (R != NULL)
            acb_one(acb_mat_entry(R, 0, 0));
        return 1;
    }

    acb_mat_init(D, n, n);
    acb_mat_init(T, n, n);
    acb_mat_init(AT, n, n);

    /* T D = A T */
    acb_mat_get_mid(T, R_approx);
    acb_mat_mul(AT, A, T, prec);

    if (acb_mat_solve(D, T, AT, prec))
    {
        acb_mat_t DD, DH;
        arb_t epsilon, eta1, eta2;

        arb_init(epsilon);
        arb_init(eta1);
        arb_init(eta2);
        acb_mat_init(DD, n, n);
        acb_mat_init(DH, n, n);

        for (i = 0; i < n; i++)
            acb_set(acb_mat_entry(DD, i, i), acb_mat_entry(D, i, i));

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                if (i != j)
                    acb_set(acb_mat_entry(DH, i, j), acb_mat_entry(D, i, j));

        diagonal_certify(epsilon, eta1, eta2, DD, DH, 2 * MAG_BITS);

        if (arb_is_finite(epsilon))
        {
            for (i = 0; i < n; i++)
            {
                /* note: paper actually uses D_c which would be better? */
                acb_set(E + i, acb_mat_entry(D, i, i));
                arb_add_error(acb_realref(E + i), eta2);
                arb_add_error(acb_imagref(E + i), eta2);
            }

            result = 1;

            /* unlikely */
            for (i = 0; i < n; i++)
                for (j = i + 1; j < n; j++)
                    if (acb_overlaps(E + i, E + j))
                        result = 0;

            if (result && (R != NULL || L != NULL))
            {
                mag_t ep, em;
                mag_init(ep);
                mag_init(em);

                arb_get_mag(ep, epsilon);
                acb_mat_zero(D);
                acb_mat_add_error_mag(D, ep);
                acb_mat_mul(D, T, D, MAG_BITS);

                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        acb_get_mag(ep, acb_mat_entry(D, i, j));
                        acb_add_error_mag(acb_mat_entry(T, i, j), ep);
                    }
                }

                if (R != NULL)
                    acb_mat_set(R, T);

                if (L != NULL)
                {
                    if (!acb_mat_inv(L, T, prec))
                        acb_mat_indeterminate(L);
                }

                mag_clear(ep);
                mag_clear(em);
            }
        }

        acb_mat_clear(DD);
        acb_mat_clear(DH);
        arb_clear(epsilon);
        arb_clear(eta1);
        arb_clear(eta2);
    }

    if (!result)
    {
        for (i = 0; i < n; i++)
            acb_indeterminate(E + i);
        if (L != NULL)
            acb_mat_indeterminate(L);
        if (R != NULL)
            acb_mat_indeterminate(R);
    }

    acb_mat_clear(D);
    acb_mat_clear(T);
    acb_mat_clear(AT);

    return result;
}
