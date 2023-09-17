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

/* Ordering is:

   [0, 'Co16'], [1, 'Co20'], [2, 'Co24'], [3, 'Co28'], [4, 'Co32'], [5,
   'Co36'], [6, 'Co38'], [7, 'Co312'], [8, 'Co40'], [9, 'Co44'], [10, 'Co46'],
   [11, 'Co410'], [12, 'Co52'], [13, 'Co54'], [14, 'Co58'], [15, 'Co60'], [16,
   'Co661'], [17, 'Co662'], [18, 'Co72'], [19, 'Co74'], [20, 'Co82'], [21,
   'Co94'], [22, 'Co100'], [23, 'Co102'], [24, 'Co122'], [25, 'Co150'] */

static void
acb_theta_g2_basic_transvectants(acb_poly_struct* res, const acb_poly_t r, slong prec)
{
    acb_poly_t s;

    acb_poly_init(s);

    /* Each block is a weight 1, 2, ..., 10, 12, 15 */
    acb_poly_set(&res[0], r);

    acb_theta_g2_transvectant(&res[1], r, r, 6, 6, 6, prec);
    acb_theta_g2_transvectant(&res[2], r, r, 6, 6, 4, prec);
    acb_theta_g2_transvectant(&res[3], r, r, 6, 6, 2, prec);

    acb_theta_g2_transvectant(&res[4], r, &res[2], 6, 4, 4, prec);
    acb_theta_g2_transvectant(&res[5], r, &res[2], 6, 4, 2, prec);
    acb_theta_g2_transvectant(&res[6], r, &res[2], 6, 4, 1, prec);
    acb_theta_g2_transvectant(&res[7], r, &res[3], 6, 8, 1, prec);

    acb_theta_g2_transvectant(&res[8], &res[2], &res[2], 4, 4, 4, prec);
    acb_theta_g2_transvectant(&res[9], r, &res[4], 6, 2, 2, prec);
    acb_theta_g2_transvectant(&res[10], r, &res[4], 6, 2, 1, prec);
    acb_theta_g2_transvectant(&res[11], &res[3], &res[2], 8, 4, 1, prec);

    acb_theta_g2_transvectant(&res[12], &res[2], &res[4], 4, 2, 2, prec);
    acb_theta_g2_transvectant(&res[13], &res[2], &res[4], 4, 2, 1, prec);
    acb_theta_g2_transvectant(&res[14], &res[3], &res[4], 8, 2, 1, prec);

    acb_theta_g2_transvectant(&res[15], &res[4], &res[4], 2, 2, 2, prec);
    acb_theta_g2_transvectant(&res[16], &res[5], &res[4], 6, 2, 1, prec);
    acb_theta_g2_transvectant(&res[17], &res[6], &res[4], 8, 2, 2, prec);

    acb_poly_mul(s, &res[4], &res[4], prec); /* C_32^2 */
    acb_theta_g2_transvectant(&res[18], r, s, 6, 4, 4, prec);
    acb_theta_g2_transvectant(&res[19], r, s, 6, 4, 3, prec);

    acb_theta_g2_transvectant(&res[20], &res[2], s, 4, 4, 3, prec);

    acb_theta_g2_transvectant(&res[21], &res[6], s, 8, 4, 4, prec);

    acb_poly_mul(s, s, &res[4], prec); /* now C_32^3 */
    acb_theta_g2_transvectant(&res[22], r, s, 6, 6, 6, prec);
    acb_theta_g2_transvectant(&res[23], r, s, 6, 6, 5, prec);

    acb_theta_g2_transvectant(&res[24], &res[6], s, 8, 6, 6, prec);

    acb_poly_mul(s, s, &res[4], prec); /* now C_32^4 */
    acb_theta_g2_transvectant(&res[25], &res[6], s, 8, 8, 8, prec);

    acb_poly_clear(s);
}

/* Ordering is:
   [0, 'Co16'], [1, 'Co20'], [2, 'Co24], [3, 'Co28'], [4, 'Co32'], [5, 'Co36'],
   [6, 'Co38'], [7, 'Co312'], [8, 'Co40'], [9, 'Co44'], [10, 'Co46'], [11,
   'Co410'], [12, 'Co52'], [13, 'Co54'], [14, 'Co58'], [15, 'Co60'], [16,
   'Co661'], [17, 'Co662'], [18, 'Co72'], [19, 'Co74'], [20, 'Co82'], [21,
   'Co94'], [22, 'Co100'], [23, 'Co102'], [24, 'Co122'], [25, 'Co150'] */

void
acb_theta_g2_basic_covariants(acb_poly_struct* res, const acb_poly_t r, slong prec)
{
    acb_t c;
    slong cofactors[ACB_THETA_G2_BASIC_NB] = {1, 60, 75, 90, 2250, 2250, 450,
        540, 11250, 67500, 13500, 13500, 168750, 67500, 405000, 10125000,
        2025000, 2700000, 151875000, 60750000, 15187500, 9112500000,
        227812500000, 13668750000, 8201250000000, 384433593750};
    slong k;

    acb_init(c);

    acb_theta_g2_basic_transvectants(res, r, prec);
    for (k = 0; k < ACB_THETA_G2_BASIC_NB; k++)
    {
        acb_set_si(c, cofactors[k]);
        acb_poly_scalar_mul(&res[k], &res[k], c, prec);
    }

    acb_clear(c);
}

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
acb_theta_g2_basic_covariants_old(acb_poly_struct* cov, const acb_poly_t r, slong prec)
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

