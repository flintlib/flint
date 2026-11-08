/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_dirichlet.h"
#include "acb_mzeta.h"

typedef int (* epstein_func)(acb_t, const acb_t, const arb_mat_t, arb_srcptr, arb_srcptr,
    int, acb_mzeta_cache_t, slong);

/* L(s) for the character mod q given by its values on 1..q-1 (q = 3, 4) */
static void
dirichlet_l_small(acb_t res, const acb_t s, ulong q, slong prec)
{
    acb_t a, t, u;
    acb_init(a);
    acb_init(t);
    acb_init(u);
    acb_set_ui(a, 1);
    acb_div_ui(a, a, q, prec);
    acb_hurwitz_zeta(t, s, a, prec);
    acb_set_ui(a, q - 1);
    acb_div_ui(a, a, q, prec);
    acb_hurwitz_zeta(u, s, a, prec);
    acb_sub(t, t, u, prec);
    acb_set_ui(a, q);
    acb_neg(u, s);
    acb_pow(a, a, u, prec);
    acb_mul(res, t, a, prec);
    acb_clear(a);
    acb_clear(t);
    acb_clear(u);
}

/* (1 - 2^(c - s)) zeta(s) */
static void
zeta_factor(acb_t res, const acb_t s, slong c, slong prec)
{
    acb_t t, u;
    acb_init(t);
    acb_init(u);
    acb_set_si(t, c);
    acb_sub(t, t, s, prec);
    acb_set_ui(u, 2);
    acb_pow(t, u, t, prec);
    acb_sub_ui(t, t, 1, prec);
    acb_neg(t, t);
    acb_zeta(u, s, prec);
    acb_mul(res, t, u, prec);
    acb_clear(t);
    acb_clear(u);
}

TEST_FUNCTION_START(acb_mzeta_epstein, state)
{
    slong iter;
    epstein_func funcs[3] = { acb_mzeta_epstein_generic, acb_mzeta_epstein_crandall,
                              acb_mzeta_epstein_crandall_dfloat };

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        slong prec, d, i, which, alg;
        arb_mat_t A;
        arb_ptr x, y;
        acb_t nu, s, ref, z, t;
        acb_mzeta_cache_t cache;
        int ok;

        which = n_randint(state, 8);
        alg = n_randint(state, 3);
        prec = (alg == 2) ? 20 + n_randint(state, 80) : 20 + n_randint(state, 200);
        d = (which <= 2) ? 1 : (which <= 5) ? 2 : 4;
        if (which == 7)
            d = 3;

        arb_mat_init(A, d, d);
        x = _arb_vec_init(d);
        y = _arb_vec_init(d);
        acb_init(nu);
        acb_init(s);
        acb_init(ref);
        acb_init(z);
        acb_init(t);
        acb_mzeta_cache_init(cache);

        arb_mat_one(A);

        /* nu: a random dyadic in [-8, 12], sometimes complex */
        arb_set_si(acb_realref(nu), (slong) n_randint(state, 5120) - 2048);
        arb_mul_2exp_si(acb_realref(nu), acb_realref(nu), -8);
        if (n_randint(state, 4) == 0)
        {
            arb_set_si(acb_imagref(nu), (slong) n_randint(state, 512) - 256);
            arb_mul_2exp_si(acb_imagref(nu), acb_imagref(nu), -6);
        }
        /* avoid the pole at nu = d of the cases with y = 0 */
        if (which == 0 || which == 1 || which == 3 || which == 5 || which == 6)
        {
            acb_sub_si(t, nu, d, prec);
            if (acb_contains_zero(t))
                acb_add_ui(nu, nu, 1, prec);
        }
        acb_mul_2exp_si(s, nu, -1);

        switch (which)
        {
            case 0:     /* 2 zeta(nu) */
                acb_zeta(ref, nu, prec + 30);
                acb_mul_2exp_si(ref, ref, 1);
                break;
            case 1:     /* x = 1/2: 2 (2^nu - 1) zeta(nu) */
                arb_set_d(x, 0.5);
                acb_set_ui(t, 2);
                acb_pow(t, t, nu, prec + 30);
                acb_sub_ui(t, t, 1, prec + 30);
                acb_zeta(ref, nu, prec + 30);
                acb_mul(ref, ref, t, prec + 30);
                acb_mul_2exp_si(ref, ref, 1);
                break;
            case 2:     /* y = 1/2: -2 eta(nu) */
                arb_set_d(y, 0.5);
                zeta_factor(ref, nu, 1, prec + 30);
                acb_mul_si(ref, ref, -2, prec + 30);
                break;
            case 3:     /* Z^2: 4 zeta(s) beta(s) */
                acb_zeta(ref, s, prec + 30);
                dirichlet_l_small(t, s, 4, prec + 30);
                acb_mul(ref, ref, t, prec + 30);
                acb_mul_2exp_si(ref, ref, 2);
                break;
            case 4:     /* Z^2, y = (1/2, 1/2): -4 eta(s) beta(s) */
                arb_set_d(y, 0.5);
                arb_set_d(y + 1, 0.5);
                zeta_factor(ref, s, 1, prec + 30);
                dirichlet_l_small(t, s, 4, prec + 30);
                acb_mul(ref, ref, t, prec + 30);
                acb_mul_si(ref, ref, -4, prec + 30);
                break;
            case 5:     /* hexagonal, Gram [2 1; 1 2]: 6 2^-s zeta(s) L(s, chi_-3) */
                arb_set_ui(arb_mat_entry(A, 0, 0), 2);
                arb_sqrt(arb_mat_entry(A, 0, 0), arb_mat_entry(A, 0, 0), prec + 100);
                arb_mul_2exp_si(arb_mat_entry(A, 0, 1), arb_mat_entry(A, 0, 0), -1);
                arb_set_ui(arb_mat_entry(A, 1, 1), 3);
                arb_mul_2exp_si(arb_mat_entry(A, 1, 1), arb_mat_entry(A, 1, 1), -1);
                arb_sqrt(arb_mat_entry(A, 1, 1), arb_mat_entry(A, 1, 1), prec + 100);
                acb_zeta(ref, s, prec + 30);
                dirichlet_l_small(t, s, 3, prec + 30);
                acb_mul(ref, ref, t, prec + 30);
                acb_set_ui(t, 2);
                acb_neg(z, s);
                acb_pow(t, t, z, prec + 30);
                acb_mul(ref, ref, t, prec + 30);
                acb_mul_ui(ref, ref, 6, prec + 30);
                break;
            case 6:     /* Z^4: 8 (1 - 4^(1-s)) zeta(s) zeta(s-1) */
                acb_set_si(t, 2);
                acb_mul_2exp_si(z, s, 1);
                acb_sub(t, t, z, prec + 30);
                acb_set_ui(z, 2);
                acb_pow(t, z, t, prec + 30);
                acb_sub_ui(t, t, 1, prec + 30);
                acb_neg(t, t);
                acb_zeta(ref, s, prec + 30);
                acb_mul(ref, ref, t, prec + 30);
                acb_sub_ui(t, s, 1, prec + 30);
                acb_zeta(t, t, prec + 30);
                acb_mul(ref, ref, t, prec + 30);
                acb_mul_ui(ref, ref, 8, prec + 30);
                break;
            default:    /* Madelung: nu = 1, y = (1/2)^3 */
                acb_one(nu);
                for (i = 0; i < 3; i++)
                    arb_set_d(y + i, 0.5 + (slong) n_randint(state, 5) - 2);
                arb_set_str(acb_realref(ref), "-1.747564594633182190636212035544397403485161436624741758152825350765040623532761179890758362694607889930832581538753710593282029944183828013036933002156599363282376607172297568659238037167203810410603 +/- 1e-197", 700);
                arb_zero(acb_imagref(ref));
                prec = FLINT_MIN(prec, 600);
                break;
        }

        ok = funcs[alg](z, nu, A, x, y, 0, n_randint(state, 2) ? cache : NULL, prec);

        if (ok && !acb_overlaps(z, ref))
        {
            flint_printf("FAIL: overlap (case %wd, algorithm %wd)\n", which, alg);
            flint_printf("prec = %wd, d = %wd\n", prec, d);
            flint_printf("nu = "); acb_printd(nu, 20); flint_printf("\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
            flint_printf("ref = "); acb_printd(ref, 30); flint_printf("\n");
            flint_abort();
        }

        if (ok && acb_is_finite(ref) && acb_rel_accuracy_bits(ref) > prec + 10 && !acb_is_zero(ref)
            && acb_rel_accuracy_bits(z) < prec - 40 - 2 * arf_abs_bound_lt_2exp_si(arb_midref(acb_realref(nu))))
        {
            flint_printf("FAIL: accuracy (case %wd, algorithm %wd)\n", which, alg);
            flint_printf("prec = %wd, d = %wd\n", prec, d);
            flint_printf("nu = "); acb_printd(nu, 20); flint_printf("\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
            flint_printf("ref = "); acb_printd(ref, 30); flint_printf("\n");
            flint_abort();
        }

        arb_mat_clear(A);
        _arb_vec_clear(x, d);
        _arb_vec_clear(y, d);
        acb_clear(nu);
        acb_clear(s);
        acb_clear(ref);
        acb_clear(z);
        acb_clear(t);
        acb_mzeta_cache_clear(cache);
    }

    TEST_FUNCTION_END(state);
}
