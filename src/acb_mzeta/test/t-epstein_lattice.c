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
#include "acb_mzeta.h"

/* exp(2 pi i t) */
static void
acb_e(acb_t res, const arb_t t, slong prec)
{
    arb_t u;
    arb_init(u);
    arb_mul_2exp_si(u, t, 1);
    arb_sin_cos_pi(acb_imagref(res), acb_realref(res), u, prec);
    arb_clear(u);
}

static void
randtest_dyadic(arb_t x, flint_rand_t state, slong bits, slong e)
{
    arb_set_si(x, (slong) n_randint(state, UWORD(1) << bits) - (WORD(1) << (bits - 1)));
    arb_mul_2exp_si(x, x, -e);
}

/* Consistency on random lattices: change of basis, translation, the
   functional equation, the regularised function, and the double-precision
   algorithm against the arb one. */
TEST_FUNCTION_START(acb_mzeta_epstein_lattice, state)
{
    slong iter;

    for (iter = 0; iter < 40 * flint_test_multiplier(); iter++)
    {
        /* the number of lattice points grows like prec^(d/2): lower
           precisions in higher dimension (the code paths are the same) */
        slong d = n_randint(state, 8) ? 1 + n_randint(state, 3) : 4, i, j, k;
        slong prec = 20 + n_randint(state, (d <= 2) ? 120 : (d == 3) ? 90 : 50);
        arb_mat_t A, B, U, Ad;
        arb_ptr x, y, x2, mx;
        acb_t nu, nu2, z1, z2, t, u;
        arb_t a;
        acb_mzeta_cache_t cache;
        int test = n_randint(state, 5), ok1, ok2;

        arb_mat_init(A, d, d);
        arb_mat_init(B, d, d);
        arb_mat_init(U, d, d);
        arb_mat_init(Ad, d, d);
        x = _arb_vec_init(d);
        y = _arb_vec_init(d);
        x2 = _arb_vec_init(d);
        mx = _arb_vec_init(d);
        acb_init(nu);
        acb_init(nu2);
        acb_init(z1);
        acb_init(z2);
        acb_init(t);
        acb_init(u);
        arb_init(a);
        acb_mzeta_cache_init(cache);

        /* a random basis near the identity, exact dyadic entries */
        for (i = 0; i < d; i++)
            for (j = 0; j < d; j++)
            {
                randtest_dyadic(arb_mat_entry(A, i, j), state, 6, 7);
                if (i == j)
                    arb_add_ui(arb_mat_entry(A, i, j), arb_mat_entry(A, i, j), 1, prec);
            }

        for (i = 0; i < d; i++)
        {
            if (n_randint(state, 4))
                randtest_dyadic(x + i, state, 10, 9);
            if (n_randint(state, 4))
                randtest_dyadic(y + i, state, 10, 9);
        }

        randtest_dyadic(acb_realref(nu), state, 12, 8);
        /* complex nu (complex arithmetic for every incomplete gamma
           value) in dimension at most 3 */
        if (test != 4 && d <= 3 && n_randint(state, 3) == 0)
        {
            randtest_dyadic(acb_imagref(nu), state, 8, 6);
            if (d == 3)
                prec = FLINT_MIN(prec, 60);
        }

        if (test == 0)
        {
            /* unimodular change of basis A -> A U */
            arb_mat_one(U);
            for (k = 0; k < 3 * d; k++)
            {
                slong r = n_randint(state, d), c = n_randint(state, d), m = (slong) n_randint(state, 5) - 2;
                if (r != c)
                    for (i = 0; i < d; i++)
                        arb_addmul_si(arb_mat_entry(U, i, c), arb_mat_entry(U, i, r), m, prec);
            }
            arb_mat_mul(B, A, U, 10 * prec);
            ok1 = acb_mzeta_epstein_crandall(z1, nu, A, x, y, 0, cache, prec);
            ok2 = acb_mzeta_epstein_generic(z2, nu, B, x, y, 0, cache, prec);
        }
        else if (test == 1)
        {
            /* translation: Z(x + A m; y) = e(-y.A m) Z(x; y) */
            arb_zero(a);
            for (j = 0; j < d; j++)
                arb_set_si(arb_mat_entry(U, j, 0), (slong) n_randint(state, 7) - 3);
            for (i = 0; i < d; i++)
            {
                arb_set(x2 + i, x + i);
                for (j = 0; j < d; j++)
                    arb_addmul(x2 + i, arb_mat_entry(A, i, j), arb_mat_entry(U, j, 0), 10 * prec);
            }
            for (i = 0; i < d; i++)
                arb_submul(a, y + i, x2 + i, 10 * prec), arb_addmul(a, y + i, x + i, 10 * prec);
            ok1 = acb_mzeta_epstein_crandall(z1, nu, A, x, y, 0, cache, prec);
            ok2 = acb_mzeta_epstein(z2, nu, A, x2, y, 0, prec);   /* thread-local cache */
            acb_e(t, a, 2 * prec);
            acb_mul(z1, z1, t, 2 * prec);
        }
        else if (test == 2)
        {
            /* functional equation: (V^(2/d)/pi)^(nu/2) / Gamma((d-nu)/2)
               e(x.y/2) Z_{L,nu}(x; y) is invariant under
               (L, nu, x, y) -> (L*, d - nu, y, -x) */
            arb_mat_inv(Ad, A, 4 * prec);
            arb_mat_transpose(Ad, Ad);
            for (i = 0; i < d; i++)
                arb_neg(mx + i, x + i);
            acb_set_si(nu2, d);
            acb_sub(nu2, nu2, nu, prec);
            ok1 = acb_mzeta_epstein_crandall(z1, nu, A, x, y, 0, cache, prec);
            ok2 = acb_mzeta_epstein_generic(z2, nu2, Ad, y, mx, 0, cache, prec);
            for (k = 0; k < 2; k++)
            {
                acb_ptr z = (k == 0) ? z1 : z2;
                const arb_mat_struct * M = (k == 0) ? A : Ad;
                acb_srcptr n = (k == 0) ? nu : nu2;
                arb_mat_det(a, M, 4 * prec);
                arb_abs(a, a);
                arb_root_ui(a, a, d, 4 * prec);
                arb_sqr(a, a, 4 * prec);
                arb_const_pi(acb_realref(t), 4 * prec);
                arb_div(a, a, acb_realref(t), 4 * prec);
                acb_set_arb(t, a);
                acb_mul_2exp_si(u, n, -1);
                acb_pow(t, t, u, 4 * prec);
                acb_mul(z, z, t, 4 * prec);
                acb_set_si(u, d);
                acb_sub(u, u, n, 4 * prec);
                acb_mul_2exp_si(u, u, -1);
                acb_rgamma(u, u, 4 * prec);
                acb_mul(z, z, u, 4 * prec);
            }
            arb_zero(a);
            for (i = 0; i < d; i++)
                arb_addmul(a, x + i, y + i, 4 * prec);
            arb_mul_2exp_si(a, a, -1);
            acb_e(t, a, 4 * prec);
            acb_mul(z1, z1, t, 4 * prec);
            arb_neg(a, a);
            acb_e(t, a, 4 * prec);
            acb_mul(z2, z2, t, 4 * prec);
        }
        else if (test == 3)
        {
            /* regularised: e(x.y) Z - s_nu(y)/V with
               s_nu(y) = pi^(nu/2) Gamma((d-nu)/2) / Gamma(nu/2) (pi y^2)^((nu-d)/2) */
            arb_t y2;
            arb_init(y2);
            for (i = 0; i < d; i++)
                arb_addmul(y2, y + i, y + i, 4 * prec);
            ok1 = acb_mzeta_epstein_crandall(z1, nu, A, x, y, ACB_MZETA_EPSTEIN_REGULARIZED, cache, prec);
            ok2 = acb_mzeta_epstein_generic(z2, nu, A, x, y, 0, cache, prec);
            if (arb_is_positive(y2))
            {
                arb_zero(a);
                for (i = 0; i < d; i++)
                    arb_addmul(a, x + i, y + i, 4 * prec);
                acb_e(t, a, 4 * prec);
                acb_mul(z2, z2, t, 4 * prec);
                /* s_nu(y) / V */
                arb_const_pi(a, 4 * prec);
                arb_mul(y2, y2, a, 4 * prec);
                acb_set_arb(t, y2);
                acb_sub_si(u, nu, d, 4 * prec);
                acb_mul_2exp_si(u, u, -1);
                acb_pow(t, t, u, 4 * prec);
                acb_neg(u, u);
                acb_gamma(u, u, 4 * prec);
                acb_mul(t, t, u, 4 * prec);
                acb_mul_2exp_si(u, nu, -1);
                acb_rgamma(u, u, 4 * prec);
                acb_mul(t, t, u, 4 * prec);
                acb_set_arb(u, a);
                acb_mul_2exp_si(nu2, nu, -1);
                acb_pow(u, u, nu2, 4 * prec);
                acb_mul(t, t, u, 4 * prec);
                arb_mat_det(a, A, 4 * prec);
                arb_abs(a, a);
                acb_div_arb(t, t, a, 4 * prec);
                acb_sub(z2, z2, t, 4 * prec);
            }
            else
                acb_set(z2, z1);   /* y = 0: nothing to compare */
            arb_clear(y2);
        }
        else
        {
            /* double precision against arb (real nu); sometimes a nu
               with low-order bits (nu/2 and (d - nu)/2 are then not both
               doubles) */
            prec = 20 + n_randint(state, (d <= 3) ? 76 : 40);
            if (n_randint(state, 4) == 0)
            {
                arb_set_d(acb_realref(nu), (double) n_randint(state, 1000) / (1 + n_randint(state, 999)));
                arb_zero(acb_imagref(nu));
            }
            ok1 = acb_mzeta_epstein_crandall_dfloat(z1, nu, A, x, y,
                n_randint(state, 2) ? ACB_MZETA_EPSTEIN_REGULARIZED : 0, cache, prec);
            ok2 = acb_mzeta_epstein_crandall(z2, nu, A, x, y,
                0, cache, prec);
            if (ok1)
            {
                /* same flags for both */
                ok2 = acb_mzeta_epstein_crandall(z2, nu, A, x, y, 0, cache, prec);
                ok1 = acb_mzeta_epstein_crandall_dfloat(z1, nu, A, x, y, 0, cache, prec);
            }
        }

        if (ok1 && ok2 && !acb_overlaps(z1, z2))
        {
            flint_printf("FAIL (test %d)\n", test);
            flint_printf("d = %wd, prec = %wd\n", d, prec);
            flint_printf("A = "); arb_mat_printd(A, 10); flint_printf("\n");
            flint_printf("x = "); _arb_vec_printd(x, d, 10); flint_printf("\n");
            flint_printf("y = "); _arb_vec_printd(y, d, 10); flint_printf("\n");
            flint_printf("nu = "); acb_printd(nu, 10); flint_printf("\n");
            flint_printf("z1 = "); acb_printd(z1, 30); flint_printf("\n");
            flint_printf("z2 = "); acb_printd(z2, 30); flint_printf("\n");
            flint_abort();
        }

        arb_mat_clear(A);
        arb_mat_clear(B);
        arb_mat_clear(U);
        arb_mat_clear(Ad);
        _arb_vec_clear(x, d);
        _arb_vec_clear(y, d);
        _arb_vec_clear(x2, d);
        _arb_vec_clear(mx, d);
        acb_clear(nu);
        acb_clear(nu2);
        acb_clear(z1);
        acb_clear(z2);
        acb_clear(t);
        acb_clear(u);
        arb_clear(a);
        acb_mzeta_cache_clear(cache);
    }

    TEST_FUNCTION_END(state);
}
