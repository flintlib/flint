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
arb_mat_inf_norm(arb_t res, const arb_mat_t mat, slong prec)
{
    arb_t t, u;
    slong k, j;

    arb_init(t);
    arb_init(u);

    arb_zero(res);
    for (k = 0; k < arb_mat_nrows(mat); k++)
    {
        arb_zero(t);
        for (j = 0; j < arb_mat_ncols(mat); j++)
        {
            arb_abs(u, arb_mat_entry(mat, k, j));
            arb_max(t, t, u, prec);
        }
        arb_max(res, res, t, prec);
    }

    arb_clear(t);
    arb_clear(u);
}

static void
_arb_vec_inf_norm(arb_t res, arb_srcptr v, slong nb, slong prec)
{
    arb_t t;
    slong k;

    arb_init(t);

    arb_zero(res);
    for (k = 0; k < nb; k++)
    {
        arb_abs(t, &v[k]);
        arb_max(res, res, t, prec);
    }

    arb_clear(t);
}

/* Evaluate upper bound on the tail */
static void
acb_theta_jet_naive_tail(arb_t res, const arf_t R2, const arb_mat_t C, arb_srcptr v, slong ord)
{
    slong g = arb_mat_nrows(C);
    slong lp = ACB_THETA_LOW_PREC;
    arb_t t, u, R;
    slong k;

    arb_init(t);
    arb_init(u);
    arb_init(R);

    /* Ensure assumptions R2\geq 4, R2\geq 2*ord are satisfied */
    arb_set_arf(R, R2);
    arb_set_si(t, FLINT_MAX(4, 2 * ord));
    arb_max(R, R, t, lp);
    arb_sqrt(R, R, lp);

    /* Evaluate 2^(2*g+2) R^(g - 1) exp(-R^2) \prod(1 + gamma_i^{-1}) */
    arb_one(res);
    arb_mul_2exp_si(res, res, 2 * g + 2);

    arb_pow_ui(t, R, g - 1, lp);
    arb_mul(res, res, t, lp);

    arb_sqr(t, R, lp);
    arb_exp(t, t, lp);
    arb_mul(res, res, t, lp);

    for (k = 0; k < g; k++)
    {
        arb_inv(t, arb_mat_entry(C, k, k), lp);
        arb_add_si(t, t, 1, lp);
        arb_mul(res, res, t, lp);
    }

    /* Multiply by max(1, ||C||R + ||v||)^ord */
    arb_mat_inf_norm(t, C, lp);
    arb_mul(t, t, R, lp);
    _arb_vec_inf_norm(u, v, g, lp);
    arb_add(t, t, u, lp);
    arb_set_si(u, 1);
    arb_max(t, t, u, lp);
    arb_pow_ui(t, t, ord, lp);
    arb_mul(res, res, t, lp);

    arb_clear(t);
    arb_clear(u);
    arb_clear(R);
}

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("jet_naive_radius....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: tail is small */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 10);
        slong ord = n_randint(state, 10);
        slong prec = ACB_THETA_LOW_PREC;
        slong bits = n_randint(state, 5);
        slong exp = 10 + n_randint(state, 100);
        arb_mat_t C;
        arb_ptr v, w;
        arf_t R2, eps, t;
        arb_t bound;
        slong k;

        arb_mat_init(C, g, g);
        v = _arb_vec_init(g);
        w = _arb_vec_init(g);
        arf_init(R2);
        arf_init(eps);
        arf_init(t);
        arb_init(bound);

        arb_mat_randtest_cho(C, state, prec, bits);
        arb_mat_transpose(C, C);
        for (k = 0; k < g; k++)
        {
            arb_randtest_precise(&w[k], state, prec, bits);
        }
        arb_mat_vector_mul_col(v, C, w, prec);

        acb_theta_jet_naive_radius(R2, eps, C, v, ord, exp);
        acb_theta_jet_naive_tail(bound, R2, C, w, ord);
        arb_get_lbound_arf(t, bound, prec);

        if (arf_cmp(t, eps) > 0)
        {
            flint_printf("FAIL\n");
            arb_mat_printd(C, 5);
            flint_printf("exp = %wd, ord = %wd, eps, R2:\n", exp, ord);
            arf_printd(eps, 10);
            flint_printf("\n");
            arf_printd(R2, 10);
            flint_printf("\nbound:\n");
            arb_printd(bound, 10);
            flint_printf("\n");
            flint_abort();
        }

        arb_mat_clear(C);
        _arb_vec_clear(v, g);
        _arb_vec_clear(w, g);
        arf_clear(R2);
        arf_clear(eps);
        arf_clear(t);
        arb_clear(bound);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
