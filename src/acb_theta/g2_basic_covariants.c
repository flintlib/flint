/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"
#include "profiler.h"

/* Covariants are in 9 variables a0, ..., a6, x, y; to evaluate, we set y=1 and
   return polynomials in x as acb_poly's */
/* Ordering is: ['Co16', 'Co20', 'Co24', 'Co28', 'Co32', 'Co36', 'Co38',
   'Co312', 'Co40', 'Co44', 'Co46', 'Co410', 'Co52', 'Co54', 'Co58', 'Co60',
   'Co661', 'Co662', 'Co72', 'Co74', 'Co82', 'Co94', 'Co100', 'Co102', 'Co122',
   'Co150'] */

static char* g2_covariants_str[] = {
#include "acb_theta/g2_basic_covariants.in"
};

static void
g2_basic_covariant_eval(acb_poly_t r, const fmpz_mpoly_t cov,
    acb_srcptr chi, const fmpz_mpoly_ctx_t ctx, slong prec)
{
    slong d = fmpz_mpoly_degree_si(cov, 7, ctx);
    acb_ptr val;
    fmpz_mpoly_univar_t u;
    fmpz_mpoly_t coef;
    acb_t c;
    slong k;

    val = _acb_vec_init(9);
    fmpz_mpoly_univar_init(u, ctx);
    fmpz_mpoly_init(coef, ctx);
    acb_init(c);

    _acb_vec_set(val, chi, 7);
    acb_one(&val[7]);
    acb_one(&val[8]);
    fmpz_mpoly_to_univar(u, cov, 8, ctx);
    acb_poly_zero(r);

    for (k = 0; k <= d; k++)
    {
        fmpz_mpoly_univar_get_term_coeff(coef, u, k, ctx);
        acb_eval_fmpz_mpoly(c, coef, val, ctx, prec);
        acb_poly_set_coeff_acb(r, k, c);
    }

    _acb_vec_clear(val, 9);
    fmpz_mpoly_univar_clear(u, ctx);
    fmpz_mpoly_clear(coef, ctx);
    acb_clear(c);
}

void
acb_theta_g2_basic_covariants(acb_poly_struct* cov, const acb_poly_t r, slong prec)
{
    slong nb = ACB_THETA_G2_BASIC_NB;
    char* vars[9] = {"a0", "a1", "a2", "a3", "a4", "a5", "a6", "x", "y"};
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t pol;
    acb_ptr chi;
    slong k;
    timeit_t t0;

    fmpz_mpoly_ctx_init(ctx, 9, ORD_LEX);
    fmpz_mpoly_init(pol, ctx);
    chi = _acb_vec_init(7);

    for (k = 0; k <= 6; k++)
    {
        acb_poly_get_coeff_acb(&chi[k], r, 6 - k);
    }

    for (k = 0; k < nb; k++)
    {
        timeit_start(t0);
        fmpz_mpoly_set_str_pretty(pol, g2_covariants_str[k], (const char**) vars, ctx);
        timeit_stop(t0);
        flint_printf("reading: %wd ms\n", t0->cpu);

        timeit_start(t0);
        g2_basic_covariant_eval(&cov[k], pol, chi, ctx, prec);
        timeit_stop(t0);
        flint_printf("eval: %wd ms\n", t0->cpu);
    }

    fmpz_mpoly_clear(pol, ctx);
    fmpz_mpoly_ctx_clear(ctx);
    _acb_vec_clear(chi, 7);
}

