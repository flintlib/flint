/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "acb.h"
#include "acb_modular.h"
#include "gr.h"
#include "gr_special.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"
#include "gr_tower/lazy_impl.h"

#define MF_MUST(expr) do { int _st = (expr); if (_st != GR_SUCCESS) { flint_printf("FAIL: line %d status %d\n", __LINE__, _st); flint_abort(); } } while (0)

static void
_mf_check_overlap(gr_srcptr x, const acb_t y, const char * what, gr_srcptr tau, gr_ctx_t K)
{
    acb_t z;
    acb_init(z);
    if (gr_tower_lazy_get_acb(z, x, 128, K) != GR_SUCCESS || !acb_overlaps(z, y))
    {
        flint_printf("FAIL: %s\n", what);
        flint_printf("tau = "); gr_println(tau, K);
        flint_printf("value = "); gr_println(x, K);
        flint_printf("enclosure = "); acb_printn(z, 30, 0); flint_printf("\n");
        flint_printf("reference = "); acb_printn(y, 30, 0); flint_printf("\n");
        flint_abort();
    }
    acb_clear(z);
}

static void
_mf_check_equal(gr_srcptr x, gr_srcptr y, const char * what, gr_srcptr tau, gr_ctx_t K)
{
    truth_t eq;
    eq = gr_equal(x, y, K);
    if (eq != T_TRUE)
    {
        acb_t a, b;
        acb_init(a); acb_init(b);
        flint_printf("FAIL: %s (%d)\n", what, (int) eq);
        if (gr_tower_lazy_get_acb(a, x, 128, K) == GR_SUCCESS && gr_tower_lazy_get_acb(b, y, 128, K) == GR_SUCCESS)
        {
            acb_printn(a, 30, 0); flint_printf("\n");
            acb_printn(b, 30, 0); flint_printf("\n");
        }
        acb_clear(a); acb_clear(b);
        flint_printf("tau = "); gr_println(tau, K);
        flint_printf("x = "); gr_println(x, K);
        flint_printf("y = "); gr_println(y, K);
        flint_abort();
    }
}

static void
_mf_check_str(gr_srcptr x, const char * s, const char * what, gr_ctx_t K)
{
    gr_ptr y;
    GR_TMP_INIT(y, K);
    MF_MUST(gr_set_str(y, s, K));
    if (gr_equal(x, y, K) != T_TRUE)
    {
        flint_printf("FAIL: %s = %s\n", what, s);
        flint_printf("x = "); gr_println(x, K);
        flint_abort();
    }
    GR_TMP_CLEAR(y, K);
}

TEST_FUNCTION_START(gr_tower_modforms, state)
{
    gr_ctx_t QQ, K;
    slong iter;

    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);

    /* values at random points against acb, and exact identities (in a
       fresh context for each point: the points of earlier iterations
       would otherwise become anchors of the later ones, at random levels,
       which the tests of commensurable points below cover) */
    for (iter = 0; iter < 3 * flint_test_multiplier(); iter++)
    {
        gr_ptr tau, t1, t2, t3, t4, x, y, z, w;
        acb_t at, r1, r2, r3, r4, zero, ref;
        fmpq_t c;
        slong k;

        gr_ctx_clear(K);
        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);

        GR_TMP_INIT5(tau, t1, t2, t3, t4, K);
        GR_TMP_INIT4(x, y, z, w, K);
        acb_init(at); acb_init(r1); acb_init(r2); acb_init(r3); acb_init(r4); acb_init(zero); acb_init(ref);
        fmpq_init(c);

        /* tau = r + i s with r rational and s rational, a square root or pi */
        fmpz_set_si(fmpq_numref(c), (slong) n_randint(state, 21) - 10);
        fmpz_set_si(fmpq_denref(c), 1 + n_randint(state, 6));
        fmpq_canonicalise(c);
        MF_MUST(gr_set_fmpq(tau, c, K));
        k = 1 + n_randint(state, 6);
        switch (n_randint(state, 3))
        {
            case 0: MF_MUST(gr_set_si(x, k, K)); MF_MUST(gr_div_ui(x, x, 1 + n_randint(state, 4), K)); break;
            case 1: MF_MUST(gr_set_si(x, k + 1, K)); MF_MUST(gr_sqrt(x, x, K)); MF_MUST(gr_div_ui(x, x, 1 + n_randint(state, 3), K)); break;
            default: MF_MUST(gr_pi(x, K)); MF_MUST(gr_div_ui(x, x, k, K)); break;
        }
        MF_MUST(gr_i(y, K));
        MF_MUST(gr_mul(x, x, y, K));
        MF_MUST(gr_add(tau, tau, x, K));
        MF_MUST(gr_tower_lazy_get_acb(at, tau, 256, K));

        /* lambda, j, eta, theta constants, Eisenstein series */
        MF_MUST(gr_modular_lambda(x, tau, K));
        acb_modular_lambda(ref, at, 256);
        _mf_check_overlap(x, ref, "lambda", tau, K);

        MF_MUST(gr_modular_j(y, tau, K));
        acb_modular_j(ref, at, 256);
        _mf_check_overlap(y, ref, "j", tau, K);

        MF_MUST(gr_dedekind_eta(z, tau, K));
        acb_modular_eta(ref, at, 256);
        _mf_check_overlap(z, ref, "eta", tau, K);

        MF_MUST(gr_zero(w, K));
        MF_MUST(gr_jacobi_theta(t1, t2, t3, t4, w, tau, K));
        acb_modular_theta(r1, r2, r3, r4, zero, at, 256);
        _mf_check_overlap(t2, r2, "theta_2", tau, K);
        _mf_check_overlap(t3, r3, "theta_3", tau, K);
        _mf_check_overlap(t4, r4, "theta_4", tau, K);

        /* theta_3^4 = theta_2^4 + theta_4^4, lambda = theta_2^4 / theta_3^4,
           2 eta^3 = theta_2 theta_3 theta_4 */
        MF_MUST(gr_pow_ui(w, t2, 4, K));
        MF_MUST(gr_pow_ui(t1, t4, 4, K));
        MF_MUST(gr_add(w, w, t1, K));
        MF_MUST(gr_pow_ui(t1, t3, 4, K));
        _mf_check_equal(w, t1, "Jacobi's identity", tau, K);
        MF_MUST(gr_pow_ui(w, t2, 4, K));
        MF_MUST(gr_div(w, w, t1, K));
        _mf_check_equal(w, x, "lambda in theta", tau, K);
        MF_MUST(gr_mul(w, t2, t3, K));
        MF_MUST(gr_mul(w, w, t4, K));
        MF_MUST(gr_pow_ui(t1, z, 3, K));
        MF_MUST(gr_mul_ui(t1, t1, 2, K));
        _mf_check_equal(w, t1, "theta_2 theta_3 theta_4 = 2 eta^3", tau, K);

        /* E4^3 - E6^2 = 1728 eta^24, j = E4^3 / eta^24, E8 = E4^2 */
        MF_MUST(gr_eisenstein_e(t1, 4, tau, K));
        MF_MUST(gr_eisenstein_e(t2, 6, tau, K));
        MF_MUST(gr_eisenstein_e(t3, 8, tau, K));
        acb_modular_eisenstein(r1, at, 1, 256);
        {
            arb_t zz;
            arb_init(zz);
            arb_zeta_ui(zz, 4, 256);
            arb_mul_2exp_si(zz, zz, 1);
            acb_div_arb(r1, r1, zz, 256);
            arb_clear(zz);
        }
        _mf_check_overlap(t1, r1, "E4", tau, K);
        MF_MUST(gr_pow_ui(w, t1, 3, K));
        MF_MUST(gr_pow_ui(t4, t2, 2, K));
        MF_MUST(gr_sub(w, w, t4, K));
        MF_MUST(gr_pow_ui(t4, z, 24, K));
        MF_MUST(gr_mul_ui(t4, t4, 1728, K));
        _mf_check_equal(w, t4, "E4^3 - E6^2 = 1728 Delta", tau, K);
        MF_MUST(gr_pow_ui(w, t1, 3, K));
        MF_MUST(gr_div_ui(t4, t4, 1728, K));
        MF_MUST(gr_div(w, w, t4, K));
        _mf_check_equal(w, y, "j = E4^3 / Delta", tau, K);
        MF_MUST(gr_sqr(w, t1, K));
        _mf_check_equal(w, t3, "E8 = E4^2", tau, K);

        /* E_2 against acb */
        MF_MUST(gr_eisenstein_e(t1, 2, tau, K));
        {
            acb_t h;
            acb_init(h);
            acb_set_d(h, 0.5);
            acb_elliptic_zeta(r1, h, at, 256);
            acb_mul_2exp_si(r1, r1, 1);
            acb_const_pi(h, 256);
            acb_sqr(h, h, 256);
            acb_div_ui(h, h, 3, 256);
            acb_div(r1, r1, h, 256);
            acb_clear(h);
        }
        _mf_check_overlap(t1, r1, "E2", tau, K);

        /* transformations: eta(tau + 1) = exp(pi i / 12) eta(tau),
           eta(-1/tau) = sqrt(-i tau) eta(tau), lambda(-1/tau) = 1 - lambda(tau),
           E2(-1/tau) = tau^2 E2(tau) + 6 tau / (pi i) */
        MF_MUST(gr_add_ui(w, tau, 1, K));
        MF_MUST(gr_dedekind_eta(t2, w, K));
        MF_MUST(gr_set_str(t3, "exp(pi*i/12)", K));
        MF_MUST(gr_mul(t3, t3, z, K));
        _mf_check_equal(t2, t3, "eta(tau + 1)", tau, K);
        MF_MUST(gr_inv(w, tau, K));
        MF_MUST(gr_neg(w, w, K));
        MF_MUST(gr_dedekind_eta(t2, w, K));
        MF_MUST(gr_i(t3, K));
        MF_MUST(gr_mul(t3, t3, tau, K));
        MF_MUST(gr_neg(t3, t3, K));
        MF_MUST(gr_sqrt(t3, t3, K));
        MF_MUST(gr_mul(t3, t3, z, K));
        _mf_check_equal(t2, t3, "eta(-1/tau)", tau, K);
        MF_MUST(gr_modular_lambda(t2, w, K));
        MF_MUST(gr_sub_ui(t3, x, 1, K));
        MF_MUST(gr_neg(t3, t3, K));
        _mf_check_equal(t2, t3, "lambda(-1/tau)", tau, K);
        MF_MUST(gr_eisenstein_e(t2, 2, w, K));
        MF_MUST(gr_sqr(t3, tau, K));
        MF_MUST(gr_mul(t3, t3, t1, K));
        MF_MUST(gr_set_str(t4, "6/(pi*i)", K));
        MF_MUST(gr_mul(t4, t4, tau, K));
        MF_MUST(gr_add(t3, t3, t4, K));
        _mf_check_equal(t2, t3, "E2(-1/tau)", tau, K);
        MF_MUST(gr_modular_j(t2, w, K));
        _mf_check_equal(t2, y, "j(-1/tau)", tau, K);

        GR_TMP_CLEAR5(tau, t1, t2, t3, t4, K);
        GR_TMP_CLEAR4(x, y, z, w, K);
        acb_clear(at); acb_clear(r1); acb_clear(r2); acb_clear(r3); acb_clear(r4); acb_clear(zero); acb_clear(ref);
        fmpq_clear(c);
    }

    /* theta functions of z: values against acb, quasi-periodicity, Jacobi's
       relations, the modular transformation tau -> tau + 1 */
    for (iter = 0; iter < 2 * flint_test_multiplier(); iter++)
    {
        gr_ptr tau, z, w, t[4], s[4], u, v;
        acb_t az, at, r[4];
        fmpq_t c;
        slong k;

        gr_ctx_clear(K);
        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);

        GR_TMP_INIT5(tau, z, w, u, v, K);
        GR_TMP_INIT4(t[0], t[1], t[2], t[3], K);
        GR_TMP_INIT4(s[0], s[1], s[2], s[3], K);
        acb_init(az); acb_init(at);
        for (k = 0; k < 4; k++)
            acb_init(r[k]);
        fmpq_init(c);

        /* (real parts of denominator at most 2: the modular transformations
           to the fundamental domain then have small |c tau + d|^2, whose
           prime factors p bring roots of unity of order 2p into the
           multipliers; and CM points of small conductor) */
        fmpq_set_si(c, (slong) n_randint(state, 9) - 4, 1 + n_randint(state, 2));
        MF_MUST(gr_set_fmpq(tau, c, K));
        MF_MUST(gr_set_si(u, 1 + n_randint(state, 4), K));
        if (n_randint(state, 2))
            MF_MUST(gr_sqrt(u, u, K));
        MF_MUST(gr_div_ui(u, u, 1 + n_randint(state, 2), K));
        MF_MUST(gr_i(v, K));
        MF_MUST(gr_mul(u, u, v, K));
        MF_MUST(gr_add(tau, tau, u, K));

        fmpq_set_si(c, (slong) n_randint(state, 21) - 10, 1 + n_randint(state, 7));
        MF_MUST(gr_set_fmpq(z, c, K));
        fmpq_set_si(c, (slong) n_randint(state, 11) - 5, 1 + n_randint(state, 7));
        MF_MUST(gr_set_fmpq(u, c, K));
        MF_MUST(gr_mul(u, u, v, K));
        MF_MUST(gr_add(z, z, u, K));
        if (n_randint(state, 3) == 0)
        {
            MF_MUST(gr_set_ui(u, 2, K));
            MF_MUST(gr_sqrt(u, u, K));
            MF_MUST(gr_add(z, z, u, K));
        }

        MF_MUST(gr_tower_lazy_get_acb(az, z, 256, K));
        MF_MUST(gr_tower_lazy_get_acb(at, tau, 256, K));
        acb_modular_theta(r[0], r[1], r[2], r[3], az, at, 256);

        MF_MUST(gr_jacobi_theta(t[0], t[1], t[2], t[3], z, tau, K));
        for (k = 0; k < 4; k++)
            _mf_check_overlap(t[k], r[k], "theta(z, tau)", tau, K);

        /* theta_1(z + tau) = -exp(-pi i (tau + 2z)) theta_1(z), and the others */
        MF_MUST(gr_add(w, z, tau, K));
        MF_MUST(gr_jacobi_theta(s[0], s[1], s[2], s[3], w, tau, K));
        MF_MUST(gr_add(u, tau, z, K));
        MF_MUST(gr_add(u, u, z, K));
        MF_MUST(gr_set_str(v, "-pi*i", K));
        MF_MUST(gr_mul(u, u, v, K));
        MF_MUST(gr_exp(u, u, K));
        for (k = 0; k < 4; k++)
        {
            MF_MUST(gr_mul(v, u, t[k], K));
            if (k == 0 || k == 3)
                MF_MUST(gr_neg(v, v, K));
            _mf_check_equal(s[k], v, "theta(z + tau)", tau, K);
        }

        /* theta_2(z)^2 theta_4^2 = theta_4(z)^2 theta_2^2 - theta_1(z)^2 theta_3^2 */
        MF_MUST(gr_zero(w, K));
        MF_MUST(gr_jacobi_theta(s[0], s[1], s[2], s[3], w, tau, K));
        MF_MUST(gr_mul(u, t[1], s[3], K));
        MF_MUST(gr_sqr(u, u, K));
        MF_MUST(gr_mul(v, t[3], s[1], K));
        MF_MUST(gr_sqr(v, v, K));
        MF_MUST(gr_mul(w, t[0], s[2], K));
        MF_MUST(gr_sqr(w, w, K));
        MF_MUST(gr_sub(v, v, w, K));
        _mf_check_equal(u, v, "Jacobi's relation in z", tau, K);

        /* theta_3(z, tau + 1) = theta_4(z, tau) */
        MF_MUST(gr_add_ui(w, tau, 1, K));
        MF_MUST(gr_jacobi_theta_3(u, z, w, K));
        _mf_check_equal(u, t[3], "theta_3(z, tau + 1)", tau, K);

        GR_TMP_CLEAR5(tau, z, w, u, v, K);
        GR_TMP_CLEAR4(t[0], t[1], t[2], t[3], K);
        GR_TMP_CLEAR4(s[0], s[1], s[2], s[3], K);
        acb_clear(az); acb_clear(at);
        for (k = 0; k < 4; k++)
            acb_clear(r[k]);
        fmpq_clear(c);
    }

    /* aliasing of the output with the inputs (the parser evaluates
       f(x) in place) */
    {
        gr_ptr tau, z, x, y, w;
        int k;
        GR_TMP_INIT5(tau, z, x, y, w, K);
        MF_MUST(gr_set_str(tau, "1/3 + pi*i/4", K));
        MF_MUST(gr_set_str(z, "1/5 + i/7", K));
        for (k = 0; k < 6; k++)
        {
            MF_MUST(gr_set(y, tau, K));
            switch (k)
            {
                case 0: MF_MUST(gr_modular_lambda(x, tau, K)); MF_MUST(gr_modular_lambda(y, y, K)); break;
                case 1: MF_MUST(gr_modular_j(x, tau, K)); MF_MUST(gr_modular_j(y, y, K)); break;
                case 2: MF_MUST(gr_dedekind_eta(x, tau, K)); MF_MUST(gr_dedekind_eta(y, y, K)); break;
                case 3: MF_MUST(gr_modular_delta(x, tau, K)); MF_MUST(gr_modular_delta(y, y, K)); break;
                case 4: MF_MUST(gr_eisenstein_e(x, 2, tau, K)); MF_MUST(gr_eisenstein_e(y, 2, y, K)); break;
                default: MF_MUST(gr_eisenstein_g(x, 4, tau, K)); MF_MUST(gr_eisenstein_g(y, 4, y, K)); break;
            }
            _mf_check_equal(x, y, "aliasing (tau)", tau, K);
        }
        for (k = 0; k < 2; k++)
        {
            /* theta_1(z, tau) with res = z, then res = tau */
            MF_MUST(gr_jacobi_theta_1(x, z, tau, K));
            if (k == 0)
            {
                MF_MUST(gr_set(y, z, K));
                MF_MUST(gr_jacobi_theta_1(y, y, tau, K));
            }
            else
            {
                MF_MUST(gr_set(y, tau, K));
                MF_MUST(gr_jacobi_theta_1(y, z, y, K));
            }
            _mf_check_equal(x, y, "aliasing (theta)", tau, K);
        }
        /* jacobi_theta with res1 = z, res2 = tau */
        {
            gr_ptr t3, t4;
            GR_TMP_INIT2(t3, t4, K);
            MF_MUST(gr_set(x, z, K));
            MF_MUST(gr_set(y, tau, K));
            MF_MUST(gr_jacobi_theta(x, y, t3, t4, x, y, K));
            MF_MUST(gr_jacobi_theta_2(w, z, tau, K));
            _mf_check_equal(w, y, "aliasing (theta, four outputs)", tau, K);
            GR_TMP_CLEAR2(t3, t4, K);
        }
        /* 2F1(a, b; c; z) with res = z, res = a; K(m) with res = m */
        MF_MUST(gr_set_str(w, "1/3", K));
        MF_MUST(gr_hypgeom_2f1(x, w, w, tau, z, 0, K));
        MF_MUST(gr_set(y, z, K));
        MF_MUST(gr_hypgeom_2f1(y, w, w, tau, y, 0, K));
        _mf_check_equal(x, y, "aliasing (2F1, z)", tau, K);
        MF_MUST(gr_set(y, w, K));
        MF_MUST(gr_hypgeom_2f1(y, y, w, tau, z, 0, K));
        _mf_check_equal(x, y, "aliasing (2F1, a)", tau, K);
        MF_MUST(gr_elliptic_k(x, z, K));
        MF_MUST(gr_set(y, z, K));
        MF_MUST(gr_elliptic_k(y, y, K));
        _mf_check_equal(x, y, "aliasing (K)", tau, K);
        GR_TMP_CLEAR5(tau, z, x, y, w, K);
    }

    /* the cached values at a point survive the evaluations at other
       points which the computation of further values there involves (E4
       caches theta_3 and lambda at tau0; theta_2, theta_4 then evaluate
       lambda at the conjugate point, whose caching once evicted theta_3
       at tau0 while it was being read): eta and the theta constants
       after E4, in a fresh context */
    {
        gr_ctx_t K6;
        gr_ptr tau, x, y, t[4];
        acb_t at, r1, r2, r3, r4, zero;

        gr_ctx_init_tower_lazy(K6, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_TMP_INIT3(tau, x, y, K6);
        GR_TMP_INIT4(t[0], t[1], t[2], t[3], K6);
        acb_init(at); acb_init(r1); acb_init(r2); acb_init(r3); acb_init(r4); acb_init(zero);

        MF_MUST(gr_set_str(tau, "7/6 + 5*i/4", K6));
        MF_MUST(gr_tower_lazy_get_acb(at, tau, 128, K6));
        MF_MUST(gr_eisenstein_e(x, 4, tau, K6));
        MF_MUST(gr_dedekind_eta(y, tau, K6));
        acb_modular_eta(r1, at, 128);
        _mf_check_overlap(y, r1, "eta after E4", tau, K6);
        MF_MUST(gr_zero(x, K6));
        MF_MUST(gr_jacobi_theta(t[0], t[1], t[2], t[3], x, tau, K6));
        acb_modular_theta(r1, r2, r3, r4, zero, at, 128);
        _mf_check_overlap(t[1], r2, "theta_2 after E4", tau, K6);
        _mf_check_overlap(t[2], r3, "theta_3 after E4", tau, K6);
        _mf_check_overlap(t[3], r4, "theta_4 after E4", tau, K6);

        GR_TMP_CLEAR3(tau, x, y, K6);
        GR_TMP_CLEAR4(t[0], t[1], t[2], t[3], K6);
        acb_clear(at); acb_clear(r1); acb_clear(r2); acb_clear(r3); acb_clear(r4); acb_clear(zero);
        gr_ctx_clear(K6);
    }

    /* the stabilizer of tau0 = i, rho acting on z: theta_3(-i z, i) =
       exp(pi z^2) theta_3(z, i), theta_1(rho z, rho) against acb */
    {
        gr_ptr tau, z, w, x, y;
        acb_t az, at, r1, r2, r3, r4;
        GR_TMP_INIT5(tau, z, w, x, y, K);
        acb_init(az); acb_init(at); acb_init(r1); acb_init(r2); acb_init(r3); acb_init(r4);

        MF_MUST(gr_set_str(tau, "i", K));
        MF_MUST(gr_set_str(z, "(1 + i)/3", K));
        MF_MUST(gr_jacobi_theta_3(x, z, tau, K));
        MF_MUST(gr_set_str(w, "-i", K));
        MF_MUST(gr_mul(w, w, z, K));
        MF_MUST(gr_jacobi_theta_3(y, w, tau, K));
        MF_MUST(gr_sqr(w, z, K));
        MF_MUST(gr_set_str(tau, "pi", K));
        MF_MUST(gr_mul(w, w, tau, K));
        MF_MUST(gr_exp(w, w, K));
        MF_MUST(gr_mul(x, x, w, K));
        MF_MUST(gr_set_str(tau, "i", K));
        _mf_check_equal(x, y, "theta_3(-i z, i)", tau, K);

        MF_MUST(gr_set_str(tau, "(1+sqrt(-3))/2", K));
        MF_MUST(gr_mul(w, z, tau, K));
        MF_MUST(gr_jacobi_theta_1(x, w, tau, K));
        MF_MUST(gr_tower_lazy_get_acb(az, w, 128, K));
        MF_MUST(gr_tower_lazy_get_acb(at, tau, 128, K));
        acb_modular_theta(r1, r2, r3, r4, az, at, 128);
        _mf_check_overlap(x, r1, "theta_1(rho z, rho)", tau, K);

        GR_TMP_CLEAR5(tau, z, w, x, y, K);
        acb_clear(az); acb_clear(at); acb_clear(r1); acb_clear(r2); acb_clear(r3); acb_clear(r4);
    }

    /* commensurable points (level 2^k): the modular equations hold */
    {
        gr_ptr tau, t2, x, y, z, u, v, w, p;
        acb_t a, b;
        static const char * phi2[11][3] = {
            {"1", "3", "0"}, {"1", "0", "3"}, {"-1", "2", "2"}, {"1488", "2", "1"}, {"1488", "1", "2"},
            {"-162000", "2", "0"}, {"-162000", "0", "2"}, {"40773375", "1", "1"},
            {"8748000000", "1", "0"}, {"8748000000", "0", "1"}, {"-157464000000000", "0", "0"}};
        slong k;

        GR_TMP_INIT5(tau, t2, x, y, z, K);
        GR_TMP_INIT4(u, v, w, p, K);
        acb_init(a); acb_init(b);

        MF_MUST(gr_set_str(tau, "pi*i/3 + 1/5", K));
        MF_MUST(gr_mul_ui(t2, tau, 2, K));

        /* Landen: lambda(2 tau) = ((1 - k') / (1 + k'))^2, k' = theta_4^2 / theta_3^2 */
        MF_MUST(gr_modular_lambda(y, t2, K));
        MF_MUST(gr_zero(z, K));
        MF_MUST(gr_jacobi_theta_3(w, z, tau, K));
        MF_MUST(gr_jacobi_theta_4(u, z, tau, K));
        MF_MUST(gr_div(u, u, w, K));
        MF_MUST(gr_sqr(u, u, K));
        MF_MUST(gr_sub_ui(v, u, 1, K));
        MF_MUST(gr_neg(v, v, K));
        MF_MUST(gr_add_ui(u, u, 1, K));
        MF_MUST(gr_div(v, v, u, K));
        MF_MUST(gr_sqr(v, v, K));
        _mf_check_equal(y, v, "lambda(2 tau) (Landen)", tau, K);

        /* the modular polynomial Phi_2(j(tau), j(2 tau)) = 0 */
        MF_MUST(gr_modular_j(x, tau, K));
        MF_MUST(gr_modular_j(y, t2, K));
        MF_MUST(gr_zero(u, K));
        for (k = 0; k < 11; k++)
        {
            MF_MUST(gr_pow_ui(p, x, phi2[k][1][0] - '0', K));
            MF_MUST(gr_pow_ui(w, y, phi2[k][2][0] - '0', K));
            MF_MUST(gr_mul(p, p, w, K));
            MF_MUST(gr_set_str(w, phi2[k][0], K));
            MF_MUST(gr_mul(p, p, w, K));
            MF_MUST(gr_add(u, u, p, K));
        }
        _mf_check_str(u, "0", "Phi_2(j(tau), j(2 tau))", K);

        /* theta_4(2 tau) = eta(tau)^2 / eta(2 tau) */
        MF_MUST(gr_dedekind_eta(x, tau, K));
        MF_MUST(gr_dedekind_eta(y, t2, K));
        MF_MUST(gr_sqr(x, x, K));
        MF_MUST(gr_div(x, x, y, K));
        MF_MUST(gr_jacobi_theta_4(y, z, t2, K));
        _mf_check_equal(x, y, "theta_4(2 tau) = eta(tau)^2 / eta(2 tau)", tau, K);

        /* 2 E_2(2 tau) - E_2(tau) = theta_3(2 tau)^4 + theta_2(2 tau)^4 */
        MF_MUST(gr_eisenstein_e(x, 2, tau, K));
        MF_MUST(gr_eisenstein_e(y, 2, t2, K));
        MF_MUST(gr_mul_ui(y, y, 2, K));
        MF_MUST(gr_sub(y, y, x, K));
        MF_MUST(gr_jacobi_theta_3(u, z, t2, K));
        MF_MUST(gr_jacobi_theta_2(v, z, t2, K));
        MF_MUST(gr_pow_ui(u, u, 4, K));
        MF_MUST(gr_pow_ui(v, v, 4, K));
        MF_MUST(gr_add(u, u, v, K));
        _mf_check_equal(y, u, "E_2 duplication", tau, K);

        /* values at (3 tau + 1)/4 and -1/(4 tau) against acb */
        MF_MUST(gr_mul_ui(u, tau, 3, K));
        MF_MUST(gr_add_ui(u, u, 1, K));
        MF_MUST(gr_div_ui(u, u, 4, K));
        MF_MUST(gr_modular_j(x, u, K));
        MF_MUST(gr_tower_lazy_get_acb(a, u, 128, K));
        acb_modular_j(b, a, 128);
        _mf_check_overlap(x, b, "j((3 tau + 1)/4)", tau, K);
        MF_MUST(gr_dedekind_eta(x, u, K));
        acb_modular_eta(b, a, 128);
        _mf_check_overlap(x, b, "eta((3 tau + 1)/4)", tau, K);
        MF_MUST(gr_mul_ui(u, tau, 4, K));
        MF_MUST(gr_inv(u, u, K));
        MF_MUST(gr_neg(u, u, K));
        MF_MUST(gr_modular_lambda(x, u, K));
        MF_MUST(gr_tower_lazy_get_acb(a, u, 128, K));
        acb_modular_lambda(b, a, 128);
        _mf_check_overlap(x, b, "lambda(-1/(4 tau))", tau, K);

        GR_TMP_CLEAR5(tau, t2, x, y, z, K);
        GR_TMP_CLEAR4(u, v, w, p, K);
        acb_clear(a); acb_clear(b);
    }

    /* commensurable points (level 2^a 3^b): with the Hauptmodul
       t = (eta(tau) / eta(3 tau))^12 of Gamma_0(3),
       j(tau) = (t + 27) (t + 243)^3 / t^3 and j(3 tau) = (t + 27) (t + 3)^3 / t
       (in a field of its own: with the anchor 2 tau of the previous
       block, the values at tau and 3 tau lie in a tower of steps for
       both halvings and the tripling, and the zero tests are costly;
       with lambda((3 tau + 1)/4) first, a root commensurable with tau
       beyond the limits of the steps, from which 3 tau has a cheaper
       chain than from tau: the values at 3 tau come from the root of
       tau, the most recently used) */
    {
        gr_ctx_t K3;
        gr_ptr tau, x, y, t, u, v;
        gr_ctx_t CC;
        acb_t a, b;
        slong k;

        gr_ctx_init_tower_lazy(K3, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_TMP_INIT4(tau, x, y, t, K3);
        GR_TMP_INIT2(u, v, K3);
        acb_init(a); acb_init(b);
        gr_ctx_init_complex_acb(CC, 128);

        MF_MUST(gr_set_str(tau, "pi*i/3 + 1/5", K3));
        MF_MUST(gr_mul_ui(u, tau, 3, K3));
        MF_MUST(gr_add_ui(u, u, 1, K3));
        MF_MUST(gr_div_ui(u, u, 4, K3));
        MF_MUST(gr_modular_lambda(x, u, K3));
        MF_MUST(gr_mul_ui(u, tau, 3, K3));
        MF_MUST(gr_dedekind_eta(x, tau, K3));
        MF_MUST(gr_dedekind_eta(y, u, K3));
        MF_MUST(gr_div(t, x, y, K3));
        MF_MUST(gr_pow_ui(t, t, 12, K3));

        MF_MUST(gr_modular_j(x, tau, K3));
        MF_MUST(gr_pow_ui(y, t, 3, K3));
        MF_MUST(gr_mul(x, x, y, K3));
        MF_MUST(gr_add_ui(y, t, 243, K3));
        MF_MUST(gr_pow_ui(y, y, 3, K3));
        MF_MUST(gr_add_ui(v, t, 27, K3));
        MF_MUST(gr_mul(y, y, v, K3));
        _mf_check_equal(x, y, "j(tau) in the Hauptmodul of level 3", tau, K3);

        MF_MUST(gr_modular_j(x, u, K3));
        MF_MUST(gr_mul(x, x, t, K3));
        MF_MUST(gr_add_ui(y, t, 3, K3));
        MF_MUST(gr_pow_ui(y, y, 3, K3));
        MF_MUST(gr_mul(y, y, v, K3));
        _mf_check_equal(x, y, "j(3 tau) in the Hauptmodul of level 3", tau, K3);

        /* the same at tau / 3 (going down) */
        MF_MUST(gr_div_ui(u, tau, 3, K3));
        MF_MUST(gr_dedekind_eta(x, u, K3));
        MF_MUST(gr_dedekind_eta(y, tau, K3));
        MF_MUST(gr_div(t, x, y, K3));
        MF_MUST(gr_pow_ui(t, t, 12, K3));
        MF_MUST(gr_modular_j(x, tau, K3));
        MF_MUST(gr_mul(x, x, t, K3));
        MF_MUST(gr_add_ui(y, t, 3, K3));
        MF_MUST(gr_pow_ui(y, y, 3, K3));
        MF_MUST(gr_add_ui(v, t, 27, K3));
        MF_MUST(gr_mul(y, y, v, K3));
        _mf_check_equal(x, y, "j(tau) = j(3 (tau / 3)) in the Hauptmodul", tau, K3);

        /* values against acb: 3 tau, (2 tau + 1)/3, 6 tau, (tau + 1)/6,
           -1/(9 tau) */
        for (k = 0; k < 5; k++)
        {
            switch (k)
            {
                case 0: MF_MUST(gr_mul_ui(u, tau, 3, K3)); break;
                case 1: MF_MUST(gr_mul_ui(u, tau, 2, K3)); MF_MUST(gr_add_ui(u, u, 1, K3)); MF_MUST(gr_div_ui(u, u, 3, K3)); break;
                case 2: MF_MUST(gr_mul_ui(u, tau, 6, K3)); break;
                case 3: MF_MUST(gr_add_ui(u, tau, 1, K3)); MF_MUST(gr_div_ui(u, u, 6, K3)); break;
                default: MF_MUST(gr_mul_ui(u, tau, 9, K3)); MF_MUST(gr_inv(u, u, K3)); MF_MUST(gr_neg(u, u, K3)); break;
            }
            /* (one function per point, E_2 at 3 tau: each is a chain of
               modular equations) */
            MF_MUST(gr_tower_lazy_get_acb(a, u, 128, K3));
            if (k == 0 || k == 3)
            {
                MF_MUST(gr_modular_lambda(x, u, K3));
                acb_modular_lambda(b, a, 128);
                _mf_check_overlap(x, b, "lambda (level 3)", tau, K3);
            }
            if (k == 1 || k == 4)
            {
                MF_MUST(gr_dedekind_eta(x, u, K3));
                acb_modular_eta(b, a, 128);
                _mf_check_overlap(x, b, "eta (level 3)", tau, K3);
            }
            if (k == 0 || k == 2)
            {
                MF_MUST(gr_eisenstein_e(x, 2, u, K3));
                MF_MUST(gr_eisenstein_e(b, 2, a, CC));
                _mf_check_overlap(x, b, "E_2 (level 3)", tau, K3);
            }
        }

        GR_TMP_CLEAR4(tau, x, y, t, K3);
        GR_TMP_CLEAR2(u, v, K3);
        acb_clear(a); acb_clear(b);
        gr_ctx_clear(CC);
        gr_ctx_clear(K3);
    }

    /* linked generators: lambda((3 tau + 1)/2) through a tripling and a
       halving from tau is a new generator linked to the chain, and 3 tau
       comes from it (the most recently used root); Jacobi's modular
       equation of degree 3 between lambda(tau) and lambda(3 tau), whose
       values come from the two roots, holds through the link (and the
       theta functions of z at the linked point against acb) */
    {
        gr_ctx_t K3;
        gr_ptr tau, x, y, l1, l2, l3, q, th[4];
        acb_t a, b, r[4];
        slong k;

        gr_ctx_init_tower_lazy(K3, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_TMP_INIT4(tau, x, y, q, K3);
        GR_TMP_INIT3(l1, l2, l3, K3);
        GR_TMP_INIT4(th[0], th[1], th[2], th[3], K3);
        acb_init(a); acb_init(b);
        for (k = 0; k < 4; k++)
            acb_init(r[k]);

        MF_MUST(gr_set_str(tau, "pi*i/3 + 1/5", K3));
        MF_MUST(gr_modular_lambda(l1, tau, K3));
        MF_MUST(gr_mul_ui(x, tau, 3, K3));
        MF_MUST(gr_add_ui(x, x, 1, K3));
        MF_MUST(gr_div_ui(x, x, 2, K3));
        MF_MUST(gr_modular_lambda(l2, x, K3));
        MF_MUST(gr_tower_lazy_get_acb(a, x, 128, K3));
        acb_modular_lambda(b, a, 128);
        _mf_check_overlap(l2, b, "lambda at a linked point", tau, K3);
        MF_MUST(gr_set_str(y, "1/7 + i/5", K3));
        MF_MUST(gr_jacobi_theta(th[0], th[1], th[2], th[3], y, x, K3));
        MF_MUST(gr_tower_lazy_get_acb(b, y, 128, K3));
        acb_modular_theta(r[0], r[1], r[2], r[3], b, a, 128);
        for (k = 0; k < 4; k++)
            _mf_check_overlap(th[k], r[k], "theta at a linked point", tau, K3);

        MF_MUST(gr_mul_ui(x, tau, 3, K3));
        MF_MUST(gr_modular_lambda(l3, x, K3));
        /* (lambda(tau) lambda(3 tau))^(1/4) + ((1 - lambda(tau)) (1 - lambda(3 tau)))^(1/4) = 1 */
        MF_MUST(gr_set_si(q, 1, K3));
        MF_MUST(gr_div_ui(q, q, 4, K3));
        MF_MUST(gr_mul(x, l1, l3, K3));
        MF_MUST(gr_pow(x, x, q, K3));
        MF_MUST(gr_sub_ui(y, l1, 1, K3));
        MF_MUST(gr_sub_ui(l2, l3, 1, K3));
        MF_MUST(gr_mul(y, y, l2, K3));
        MF_MUST(gr_pow(y, y, q, K3));
        MF_MUST(gr_add(x, x, y, K3));
        _mf_check_str(x, "1", "the modular equation of degree 3 through a link", K3);

        /* the records of the linked generators go with the last elements
           which involve them (their own data, in towers with their
           generators, does not keep them) */
        {
            slong n0 = LAZY(K3)->num_rebases, n1;
            MF_MUST(gr_zero(l2, K3));
            MF_MUST(gr_zero(l3, K3));
            MF_MUST(gr_zero(x, K3));
            MF_MUST(gr_zero(y, K3));
            for (k = 0; k < 4; k++)
                MF_MUST(gr_zero(th[k], K3));
            (void) gr_tower_lazy_ctx_num_towers(K3);
            n1 = LAZY(K3)->num_rebases;
            if (n0 == 0 || n1 >= n0)
            {
                flint_printf("FAIL: records of linked generators kept (%wd -> %wd)\n", n0, n1);
                gr_tower_lazy_ctx_stats(K3);
                flint_abort();
            }
        }

        GR_TMP_CLEAR4(tau, x, y, q, K3);
        GR_TMP_CLEAR3(l1, l2, l3, K3);
        GR_TMP_CLEAR4(th[0], th[1], th[2], th[3], K3);
        acb_clear(a); acb_clear(b);
        for (k = 0; k < 4; k++)
            acb_clear(r[k]);
        gr_ctx_clear(K3);
    }

    /* theta functions of z at related points: torsion points (algebraic),
       the duplication formulas in both orders of evaluation, divisions,
       the addition formulas whatever the order, and the addition theorem
       of wp */
    {
        gr_ctx_t K5, CC;
        gr_ptr tau, z, w, x, y, u, v, t[4], c[4], g2, g3;
        acb_t a, b, r1, r2, r3, r4;
        slong k, order;

        gr_ctx_init_complex_acb(CC, 128);
        acb_init(a); acb_init(b); acb_init(r1); acb_init(r2); acb_init(r3); acb_init(r4);

        for (order = 0; order < 4; order++)
        {
            gr_ctx_init_tower_lazy(K5, QQ, GR_TOWER_MERGE_EXPRESS);
            GR_TMP_INIT5(tau, z, w, x, y, K5);
            GR_TMP_INIT4(u, v, g2, g3, K5);
            GR_TMP_INIT4(t[0], t[1], t[2], t[3], K5);
            GR_TMP_INIT4(c[0], c[1], c[2], c[3], K5);

            MF_MUST(gr_set_str(tau, "1/3 + pi*i/4", K5));
            MF_MUST(gr_set_str(z, "1/5 + i/7", K5));
            MF_MUST(gr_set_str(w, "-1/9 + i/11", K5));
            MF_MUST(gr_zero(x, K5));
            MF_MUST(gr_jacobi_theta(c[0], c[1], c[2], c[3], x, tau, K5));

            if (order == 0)
            {
                /* torsion points against acb; wp(1/3) a root of psi_3 */
                static const char * tp[5] = {"1/3", "(1 + tau)/3", "1/4", "(1 + 2*tau)/6", "1/5 + 2*tau/5"};
                for (k = 0; k < 5; k++)
                {
                    (void) tp;
                    if (k == 0) MF_MUST(gr_set_str(x, "1/3", K5));
                    else if (k == 1) { MF_MUST(gr_add_ui(x, tau, 1, K5)); MF_MUST(gr_div_ui(x, x, 3, K5)); }
                    else if (k == 2) MF_MUST(gr_set_str(x, "1/4", K5));
                    else if (k == 3) { MF_MUST(gr_mul_ui(x, tau, 2, K5)); MF_MUST(gr_add_ui(x, x, 1, K5)); MF_MUST(gr_div_ui(x, x, 6, K5)); }
                    else { MF_MUST(gr_mul_ui(x, tau, 2, K5)); MF_MUST(gr_add_ui(x, x, 1, K5)); MF_MUST(gr_div_ui(x, x, 5, K5)); }
                    MF_MUST(gr_jacobi_theta(t[0], t[1], t[2], t[3], x, tau, K5));
                    MF_MUST(gr_tower_lazy_get_acb(a, x, 128, K5));
                    MF_MUST(gr_tower_lazy_get_acb(b, tau, 128, K5));
                    acb_modular_theta(r1, r2, r3, r4, a, b, 128);
                    _mf_check_overlap(t[0], r1, "theta_1 at a torsion point", tau, K5);
                    _mf_check_overlap(t[1], r2, "theta_2 at a torsion point", tau, K5);
                    _mf_check_overlap(t[2], r3, "theta_3 at a torsion point", tau, K5);
                    _mf_check_overlap(t[3], r4, "theta_4 at a torsion point", tau, K5);
                }
                MF_MUST(gr_set_str(x, "1/3", K5));
                MF_MUST(gr_weierstrass_p(u, x, tau, K5));
                MF_MUST(gr_elliptic_invariants(g2, g3, tau, K5));
                /* 3 x^4 - 3/2 g2 x^2 - 3 g3 x - g2^2/16 */
                MF_MUST(gr_pow_ui(v, u, 4, K5));
                MF_MUST(gr_mul_ui(v, v, 3, K5));
                MF_MUST(gr_sqr(y, u, K5));
                MF_MUST(gr_mul(y, y, g2, K5));
                MF_MUST(gr_mul_ui(y, y, 3, K5));
                MF_MUST(gr_div_ui(y, y, 2, K5));
                MF_MUST(gr_sub(v, v, y, K5));
                MF_MUST(gr_mul(y, g3, u, K5));
                MF_MUST(gr_mul_ui(y, y, 3, K5));
                MF_MUST(gr_sub(v, v, y, K5));
                MF_MUST(gr_sqr(y, g2, K5));
                MF_MUST(gr_div_ui(y, y, 16, K5));
                MF_MUST(gr_sub(v, v, y, K5));
                _mf_check_str(v, "0", "psi_3(wp(1/3)) = 0", K5);
            }

            /* the duplication theta_1(2z) theta_2 theta_3 theta_4 = 2 theta_1(z) ... theta_4(z),
               theta(2z) first (order 0) or theta(z) first (order 1) */
            if (order <= 1)
            {
                MF_MUST(gr_mul_ui(y, z, 2, K5));
                if (order == 0)
                {
                    MF_MUST(gr_jacobi_theta_1(x, y, tau, K5));
                    MF_MUST(gr_jacobi_theta(t[0], t[1], t[2], t[3], z, tau, K5));
                }
                else
                {
                    MF_MUST(gr_jacobi_theta(t[0], t[1], t[2], t[3], z, tau, K5));
                    MF_MUST(gr_jacobi_theta_1(x, y, tau, K5));
                }
                MF_MUST(gr_mul(x, x, c[1], K5));
                MF_MUST(gr_mul(x, x, c[2], K5));
                MF_MUST(gr_mul(x, x, c[3], K5));
                MF_MUST(gr_mul(u, t[0], t[1], K5));
                MF_MUST(gr_mul(u, u, t[2], K5));
                MF_MUST(gr_mul(u, u, t[3], K5));
                MF_MUST(gr_mul_ui(u, u, 2, K5));
                _mf_check_equal(x, u, "theta duplication", tau, K5);
            }

            /* the triplication theta_1(3z) theta_1(z) theta_4^2 =
               theta_1(2z)^2 theta_4(z)^2 - theta_4(2z)^2 theta_1(z)^2 with
               the values at 3z, 2z, z in this order (order 2: 2z and z
               divisions of 3z, made multiples of z by the rebase of 3z)
               or at 2z, z, 3z (order 3); the values at z against acb */
            if (order >= 2)
            {
                MF_MUST(gr_mul_ui(x, z, 3, K5));
                MF_MUST(gr_mul_ui(y, z, 2, K5));
                if (order == 2)
                {
                    MF_MUST(gr_jacobi_theta_1(u, x, tau, K5));
                    MF_MUST(gr_jacobi_theta_1(v, y, tau, K5));
                    MF_MUST(gr_jacobi_theta_4(g3, y, tau, K5));
                    MF_MUST(gr_jacobi_theta(t[0], t[1], t[2], t[3], z, tau, K5));
                }
                else
                {
                    MF_MUST(gr_jacobi_theta_1(v, y, tau, K5));
                    MF_MUST(gr_jacobi_theta_4(g3, y, tau, K5));
                    MF_MUST(gr_jacobi_theta(t[0], t[1], t[2], t[3], z, tau, K5));
                    MF_MUST(gr_jacobi_theta_1(u, x, tau, K5));
                }
                MF_MUST(gr_tower_lazy_get_acb(a, z, 128, K5));
                MF_MUST(gr_tower_lazy_get_acb(b, tau, 128, K5));
                acb_modular_theta(r1, r2, r3, r4, a, b, 128);
                _mf_check_overlap(t[0], r1, "theta_1 after the rebase", tau, K5);
                _mf_check_overlap(t[1], r2, "theta_2 after the rebase", tau, K5);
                _mf_check_overlap(t[2], r3, "theta_3 after the rebase", tau, K5);
                _mf_check_overlap(t[3], r4, "theta_4 after the rebase", tau, K5);
                MF_MUST(gr_mul(u, u, t[0], K5));
                MF_MUST(gr_mul(u, u, c[3], K5));
                MF_MUST(gr_mul(u, u, c[3], K5));
                MF_MUST(gr_mul(v, v, t[3], K5));
                MF_MUST(gr_sqr(v, v, K5));
                MF_MUST(gr_mul(g3, g3, t[0], K5));
                MF_MUST(gr_sqr(g3, g3, K5));
                MF_MUST(gr_sub(v, v, g3, K5));
                _mf_check_equal(u, v, "theta triplication", tau, K5);
            }

            /* the addition formulas theta_k(z + w) theta_k(z - w) theta_4^2 =
               theta_k(z)^2 theta_4(w)^2 - theta_(5-k)(z)^2 theta_1(w)^2, with
               the values at z, z + w, z - w, w in this order (order 0),
               at z, w first (order 1), at w, z + w, z - w, z (order 2), or
               at z + w, z - w first (order 3: z the half of their sum, w
               of their difference); k = 4, and k = 2 in order 3 */
            MF_MUST(gr_add(x, z, w, K5));
            MF_MUST(gr_sub(y, z, w, K5));
            for (k = 4; k >= 2; k -= 2)
            {
                if (k == 2 && order != 3)
                    break;
                if (order == 0)
                {
                    MF_MUST(gr_jacobi_theta(t[0], t[1], t[2], t[3], z, tau, K5));
                    MF_MUST(gr_jacobi_theta_4(u, x, tau, K5));
                    MF_MUST(gr_jacobi_theta_4(v, y, tau, K5));
                    MF_MUST(gr_jacobi_theta_1(g2, w, tau, K5));
                    MF_MUST(gr_jacobi_theta_4(g3, w, tau, K5));
                }
                else if (order == 1)
                {
                    MF_MUST(gr_jacobi_theta(t[0], t[1], t[2], t[3], z, tau, K5));
                    MF_MUST(gr_jacobi_theta_1(g2, w, tau, K5));
                    MF_MUST(gr_jacobi_theta_4(g3, w, tau, K5));
                    MF_MUST(gr_jacobi_theta_4(u, x, tau, K5));
                    MF_MUST(gr_jacobi_theta_4(v, y, tau, K5));
                }
                else if (order == 2)
                {
                    MF_MUST(gr_jacobi_theta_1(g2, w, tau, K5));
                    MF_MUST(gr_jacobi_theta_4(u, x, tau, K5));
                    MF_MUST(gr_jacobi_theta_4(v, y, tau, K5));
                    MF_MUST(gr_jacobi_theta_4(g3, w, tau, K5));
                    MF_MUST(gr_jacobi_theta(t[0], t[1], t[2], t[3], z, tau, K5));
                }
                else
                {
                    MF_MUST((k == 4 ? gr_jacobi_theta_4 : gr_jacobi_theta_2)(u, x, tau, K5));
                    MF_MUST((k == 4 ? gr_jacobi_theta_4 : gr_jacobi_theta_2)(v, y, tau, K5));
                    MF_MUST(gr_jacobi_theta(t[0], t[1], t[2], t[3], z, tau, K5));
                    MF_MUST(gr_jacobi_theta_1(g2, w, tau, K5));
                    MF_MUST(gr_jacobi_theta_4(g3, w, tau, K5));
                }
                MF_MUST(gr_mul(u, u, v, K5));
                MF_MUST(gr_mul(u, u, c[3], K5));
                MF_MUST(gr_mul(u, u, c[3], K5));
                MF_MUST(gr_mul(v, t[k - 1], g3, K5));
                MF_MUST(gr_sqr(v, v, K5));
                MF_MUST(gr_mul(g2, t[4 - k], g2, K5));
                MF_MUST(gr_sqr(g2, g2, K5));
                MF_MUST(gr_sub(v, v, g2, K5));
                _mf_check_equal(u, v, "theta addition formula", tau, K5);
            }

            GR_TMP_CLEAR5(tau, z, w, x, y, K5);
            GR_TMP_CLEAR4(u, v, g2, g3, K5);
            GR_TMP_CLEAR4(t[0], t[1], t[2], t[3], K5);
            GR_TMP_CLEAR4(c[0], c[1], c[2], c[3], K5);
            gr_ctx_clear(K5);
        }

        /* the addition theorem of wp, wp(z + w) first */
        {
            gr_ptr P, Pz, Pzp, Pw, Pwp;
            gr_ctx_init_tower_lazy(K5, QQ, GR_TOWER_MERGE_EXPRESS);
            GR_TMP_INIT5(tau, z, w, x, y, K5);
            GR_TMP_INIT5(P, Pz, Pzp, Pw, Pwp, K5);
            MF_MUST(gr_set_str(tau, "1/3 + pi*i/4", K5));
            MF_MUST(gr_set_str(z, "1/5 + i/7", K5));
            MF_MUST(gr_set_str(w, "-1/9 + i/11", K5));
            MF_MUST(gr_add(x, z, w, K5));
            MF_MUST(gr_weierstrass_p(P, x, tau, K5));
            MF_MUST(gr_weierstrass_p(Pz, z, tau, K5));
            MF_MUST(gr_weierstrass_p_prime(Pzp, z, tau, K5));
            MF_MUST(gr_weierstrass_p(Pw, w, tau, K5));
            MF_MUST(gr_weierstrass_p_prime(Pwp, w, tau, K5));
            /* (wp'(z) - wp'(w))^2 / (4 (wp(z) - wp(w))^2) - wp(z) - wp(w) */
            MF_MUST(gr_sub(x, Pzp, Pwp, K5));
            MF_MUST(gr_sub(y, Pz, Pw, K5));
            MF_MUST(gr_div(x, x, y, K5));
            MF_MUST(gr_sqr(x, x, K5));
            MF_MUST(gr_div_ui(x, x, 4, K5));
            MF_MUST(gr_sub(x, x, Pz, K5));
            MF_MUST(gr_sub(x, x, Pw, K5));
            _mf_check_equal(P, x, "the addition theorem of wp", tau, K5);
            GR_TMP_CLEAR5(tau, z, w, x, y, K5);
            GR_TMP_CLEAR5(P, Pz, Pzp, Pw, Pwp, K5);
            gr_ctx_clear(K5);
        }

        acb_clear(a); acb_clear(b); acb_clear(r1); acb_clear(r2); acb_clear(r3); acb_clear(r4);
        gr_ctx_clear(CC);
    }

    /* Weierstrass functions: values against acb, the differential
       equation, the e_k and g_k, periodicity, the half periods */
    {
        gr_ctx_t K4, CC;
        gr_ptr tau, z, P, Pp, S, e1, e2, e3, g2, g3, u, v;
        acb_t a, b, c;

        gr_ctx_init_tower_lazy(K4, QQ, GR_TOWER_MERGE_EXPRESS);
        gr_ctx_init_complex_acb(CC, 128);
        GR_TMP_INIT4(tau, z, P, Pp, K4);
        GR_TMP_INIT4(S, e1, e2, e3, K4);
        GR_TMP_INIT4(g2, g3, u, v, K4);
        acb_init(a); acb_init(b); acb_init(c);

        MF_MUST(gr_set_str(tau, "3/10 + 11*i*sqrt(2)/10", K4));
        MF_MUST(gr_set_str(z, "21/100 + 13*i/100", K4));
        MF_MUST(gr_weierstrass_p(P, z, tau, K4));
        MF_MUST(gr_weierstrass_p_prime(Pp, z, tau, K4));
        MF_MUST(gr_weierstrass_sigma(S, z, tau, K4));
        MF_MUST(gr_elliptic_roots(e1, e2, e3, tau, K4));
        MF_MUST(gr_elliptic_invariants(g2, g3, tau, K4));

        MF_MUST(gr_tower_lazy_get_acb(a, z, 128, K4));
        MF_MUST(gr_tower_lazy_get_acb(b, tau, 128, K4));
        MF_MUST(gr_weierstrass_p(c, a, b, CC));
        _mf_check_overlap(P, c, "wp", tau, K4);
        MF_MUST(gr_weierstrass_p_prime(c, a, b, CC));
        _mf_check_overlap(Pp, c, "wp'", tau, K4);
        MF_MUST(gr_weierstrass_sigma(c, a, b, CC));
        _mf_check_overlap(S, c, "sigma", tau, K4);
        {
            acb_t d, f;
            acb_init(d); acb_init(f);
            MF_MUST(gr_elliptic_roots(c, d, f, b, CC));
            _mf_check_overlap(e1, c, "e_1", tau, K4);
            _mf_check_overlap(e2, d, "e_2", tau, K4);
            _mf_check_overlap(e3, f, "e_3", tau, K4);
            MF_MUST(gr_elliptic_invariants(c, d, b, CC));
            _mf_check_overlap(g2, c, "g_2", tau, K4);
            _mf_check_overlap(g3, d, "g_3", tau, K4);
            acb_clear(d); acb_clear(f);
        }

        /* wp'^2 = 4 wp^3 - g_2 wp - g_3 */
        MF_MUST(gr_sqr(u, Pp, K4));
        MF_MUST(gr_pow_ui(v, P, 3, K4));
        MF_MUST(gr_mul_ui(v, v, 4, K4));
        MF_MUST(gr_submul(v, g2, P, K4));
        MF_MUST(gr_sub(v, v, g3, K4));
        _mf_check_equal(u, v, "wp'^2 = 4 wp^3 - g_2 wp - g_3", tau, K4);

        /* wp'^2 = 4 (wp - e_1)(wp - e_2)(wp - e_3), g_3 = 4 e_1 e_2 e_3 */
        MF_MUST(gr_sub(v, P, e1, K4));
        MF_MUST(gr_sub(S, P, e2, K4));
        MF_MUST(gr_mul(v, v, S, K4));
        MF_MUST(gr_sub(S, P, e3, K4));
        MF_MUST(gr_mul(v, v, S, K4));
        MF_MUST(gr_mul_ui(v, v, 4, K4));
        _mf_check_equal(u, v, "wp'^2 in the e_k", tau, K4);
        MF_MUST(gr_mul(u, e1, e2, K4));
        MF_MUST(gr_mul(u, u, e3, K4));
        MF_MUST(gr_mul_ui(u, u, 4, K4));
        _mf_check_equal(u, g3, "g_3 = 4 e_1 e_2 e_3", tau, K4);

        /* periodicity and parity, wp(1/2) = e_1, wp(tau/2) = e_3 */
        MF_MUST(gr_add(u, z, tau, K4));
        MF_MUST(gr_add_ui(u, u, 1, K4));
        MF_MUST(gr_neg(u, u, K4));
        MF_MUST(gr_weierstrass_p(v, u, tau, K4));
        _mf_check_equal(v, P, "wp(-z - 1 - tau)", tau, K4);
        MF_MUST(gr_set_str(u, "1/2", K4));
        MF_MUST(gr_weierstrass_p(v, u, tau, K4));
        _mf_check_equal(v, e1, "wp(1/2) = e_1", tau, K4);
        MF_MUST(gr_div_ui(u, tau, 2, K4));
        MF_MUST(gr_weierstrass_p(v, u, tau, K4));
        _mf_check_equal(v, e3, "wp(tau/2) = e_3", tau, K4);

        /* the pole */
        MF_MUST(gr_one(u, K4));
        if (gr_weierstrass_p(v, u, tau, K4) != GR_DOMAIN)
        {
            flint_printf("FAIL: wp(1)\n");
            flint_abort();
        }

        GR_TMP_CLEAR4(tau, z, P, Pp, K4);
        GR_TMP_CLEAR4(S, e1, e2, e3, K4);
        GR_TMP_CLEAR4(g2, g3, u, v, K4);
        acb_clear(a); acb_clear(b); acb_clear(c);
        gr_ctx_clear(CC);
        gr_ctx_clear(K4);
    }

    /* CM points: algebraic values and Chowla-Selberg */
    {
        gr_ptr tau, x;
        GR_TMP_INIT2(tau, x, K);

        MF_MUST(gr_set_str(tau, "i", K));
        MF_MUST(gr_modular_j(x, tau, K));
        _mf_check_str(x, "1728", "j(i)", K);
        MF_MUST(gr_modular_lambda(x, tau, K));
        _mf_check_str(x, "1/2", "lambda(i)", K);
        MF_MUST(gr_dedekind_eta(x, tau, K));
        _mf_check_str(x, "gamma(1/4)/(2*pi^(3/4))", "eta(i)", K);
        MF_MUST(gr_eisenstein_e(x, 2, tau, K));
        _mf_check_str(x, "3/pi", "E2(i)", K);

        MF_MUST(gr_set_str(tau, "(1+sqrt(-3))/2", K));
        MF_MUST(gr_modular_j(x, tau, K));
        _mf_check_str(x, "0", "j(rho)", K);
        MF_MUST(gr_eisenstein_e(x, 4, tau, K));
        _mf_check_str(x, "0", "E4(rho)", K);

        MF_MUST(gr_set_str(tau, "sqrt(-2)", K));
        MF_MUST(gr_modular_j(x, tau, K));
        _mf_check_str(x, "8000", "j(sqrt(-2))", K);

        MF_MUST(gr_set_str(tau, "(1+sqrt(-163))/2", K));
        MF_MUST(gr_modular_j(x, tau, K));
        _mf_check_str(x, "-640320^3", "j((1+sqrt(-163))/2)", K);

        MF_MUST(gr_set_str(tau, "(5+sqrt(-163))/2", K));
        MF_MUST(gr_modular_j(x, tau, K));
        _mf_check_str(x, "-640320^3", "j((5+sqrt(-163))/2)", K);

        /* class number 2: j(sqrt(-5)) = 632000 + 282880 sqrt(5) */
        MF_MUST(gr_set_str(tau, "sqrt(-5)", K));
        MF_MUST(gr_modular_j(x, tau, K));
        _mf_check_str(x, "632000+282880*sqrt(5)", "j(sqrt(-5))", K);

        /* eta at a CM point of discriminant -7, in Gamma values */
        MF_MUST(gr_set_str(tau, "(1+sqrt(-7))/2", K));
        MF_MUST(gr_dedekind_eta(x, tau, K));
        {
            acb_t at, ref;
            acb_init(at); acb_init(ref);
            MF_MUST(gr_tower_lazy_get_acb(at, tau, 256, K));
            acb_modular_eta(ref, at, 256);
            _mf_check_overlap(x, ref, "eta((1+sqrt(-7))/2)", tau, K);
            acb_clear(at); acb_clear(ref);
        }

        /* eta(2i)^2 / eta(i)^2: algebraic (both CM) */
        MF_MUST(gr_set_str(tau, "2*i", K));
        MF_MUST(gr_dedekind_eta(x, tau, K));
        {
            gr_ptr y;
            GR_TMP_INIT(y, K);
            MF_MUST(gr_set_str(tau, "i", K));
            MF_MUST(gr_dedekind_eta(y, tau, K));
            MF_MUST(gr_div(x, x, y, K));
            MF_MUST(gr_pow_ui(x, x, 8, K));
            _mf_check_str(x, "1/8", "(eta(2i)/eta(i))^8", K);
            GR_TMP_CLEAR(y, K);
        }

        /* the lower half-plane */
        MF_MUST(gr_set_str(tau, "-i", K));
        if (gr_modular_j(x, tau, K) != GR_DOMAIN)
        {
            flint_printf("FAIL: j(-i)\n");
            flint_abort();
        }

        GR_TMP_CLEAR2(tau, x, K);
    }

    /* the reduced point of z is canonical on the boundary of the box:
       z and -z with y1 = -1/4 or x1 = -1/4 reduce to the same point, and
       a point found as a multiple of a recorded one, with a reduced
       point on the boundary, gets the values at that point */
    {
        gr_ctx_t K6;
        gr_ptr tau, z, w, x, y, u;
        acb_t az, at, r[4];
        slong k, j;
        static const char * pts[3] = {"1/8 + i/4", "-1/4 - i/8", "1/7 + i/5"};

        gr_ctx_init_tower_lazy(K6, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_TMP_INIT5(tau, z, w, x, y, K6);
        GR_TMP_INIT(u, K6);
        acb_init(az); acb_init(at);
        for (k = 0; k < 4; k++)
            acb_init(r[k]);

        MF_MUST(gr_set_str(tau, "1/5 + i", K6));
        for (j = 0; j < 3; j++)
        {
            MF_MUST(gr_set_str(z, pts[j], K6));
            MF_MUST(gr_neg(w, z, K6));
            MF_MUST(gr_jacobi_theta_1(x, z, tau, K6));
            MF_MUST(gr_jacobi_theta_1(y, w, tau, K6));
            MF_MUST(gr_neg(y, y, K6));
            _mf_check_equal(x, y, "theta_1(-z) = -theta_1(z) (boundary)", tau, K6);
        }

        /* theta(u), then theta(z) with z0 = 2u, y0 = 1/4 */
        MF_MUST(gr_set_str(u, "-1/80 + i/8", K6));
        MF_MUST(gr_set_str(z, "1/8 + i/4", K6));
        for (j = 0; j < 2; j++)
        {
            gr_ptr th[4];
            GR_TMP_INIT4(th[0], th[1], th[2], th[3], K6);
            MF_MUST(gr_jacobi_theta(th[0], th[1], th[2], th[3], (j == 0) ? u : z, tau, K6));
            MF_MUST(gr_tower_lazy_get_acb(az, (j == 0) ? u : z, 128, K6));
            MF_MUST(gr_tower_lazy_get_acb(at, tau, 128, K6));
            acb_modular_theta(r[0], r[1], r[2], r[3], az, at, 128);
            for (k = 0; k < 4; k++)
                _mf_check_overlap(th[k], r[k], "theta_k(z) after theta_k(u), boundary point", tau, K6);
            GR_TMP_CLEAR4(th[0], th[1], th[2], th[3], K6);
        }

        acb_clear(az); acb_clear(at);
        for (k = 0; k < 4; k++)
            acb_clear(r[k]);
        GR_TMP_CLEAR5(tau, z, w, x, y, K6);
        GR_TMP_CLEAR(u, K6);
        gr_ctx_clear(K6);
    }

    /* the cache of the values at the points of z is bounded: a rebase of
       an addition after the values at the points it reads would be
       evicted (29 other points between) */
    {
        gr_ctx_t K7;
        gr_ptr tau, A, P, Q, y, w, x, t, u, v, c4, lhs, rhs;
        gr_ptr ty[4], tw[4];
        acb_t az, at, r[4];
        slong k;

        gr_ctx_init_tower_lazy(K7, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_TMP_INIT5(tau, A, P, Q, y, K7);
        GR_TMP_INIT4(w, x, t, u, K7);
        GR_TMP_INIT4(v, c4, lhs, rhs, K7);
        GR_TMP_INIT4(ty[0], ty[1], ty[2], ty[3], K7);
        GR_TMP_INIT4(tw[0], tw[1], tw[2], tw[3], K7);
        acb_init(az); acb_init(at);
        for (k = 0; k < 4; k++)
            acb_init(r[k]);

        MF_MUST(gr_set_str(tau, "1/3 + pi*i/4", K7));
        MF_MUST(gr_set_str(A, "1/13 + i/17", K7));
        MF_MUST(gr_set_str(P, "1/5 + i/7", K7));
        MF_MUST(gr_set_str(Q, "-1/9 + i/11", K7));
        MF_MUST(gr_add(y, P, Q, K7));
        MF_MUST(gr_div_ui(y, y, 2, K7));
        MF_MUST(gr_sub(w, P, Q, K7));
        MF_MUST(gr_div_ui(w, w, 2, K7));

        MF_MUST(gr_jacobi_theta_1(x, A, tau, K7));
        MF_MUST(gr_jacobi_theta_1(x, P, tau, K7));
        MF_MUST(gr_jacobi_theta_1(x, Q, tau, K7));
        for (k = 1; k <= 29; k++)
        {
            char s[64];
            flint_sprintf(s, "%wd/100 + i", k);
            MF_MUST(gr_set_str(t, s, K7));
            flint_sprintf(s, "1/%wd + i/%wd", k + 13, k + 31);
            MF_MUST(gr_set_str(u, s, K7));
            MF_MUST(gr_jacobi_theta_1(x, u, t, K7));
        }

        MF_MUST(gr_jacobi_theta(ty[0], ty[1], ty[2], ty[3], y, tau, K7));
        MF_MUST(gr_tower_lazy_get_acb(az, y, 128, K7));
        MF_MUST(gr_tower_lazy_get_acb(at, tau, 128, K7));
        acb_modular_theta(r[0], r[1], r[2], r[3], az, at, 128);
        for (k = 0; k < 4; k++)
            _mf_check_overlap(ty[k], r[k], "theta_k((P + Q)/2) after evictions", tau, K7);

        /* theta_1(P) theta_1(Q) theta_4^2 = theta_1(y)^2 theta_4(w)^2 - theta_4(y)^2 theta_1(w)^2 */
        MF_MUST(gr_zero(t, K7));
        MF_MUST(gr_jacobi_theta_4(c4, t, tau, K7));
        MF_MUST(gr_jacobi_theta(tw[0], tw[1], tw[2], tw[3], w, tau, K7));
        MF_MUST(gr_jacobi_theta_1(u, P, tau, K7));
        MF_MUST(gr_jacobi_theta_1(v, Q, tau, K7));
        MF_MUST(gr_mul(lhs, u, v, K7));
        MF_MUST(gr_mul(lhs, lhs, c4, K7));
        MF_MUST(gr_mul(lhs, lhs, c4, K7));
        MF_MUST(gr_mul(rhs, ty[0], tw[3], K7));
        MF_MUST(gr_sqr(rhs, rhs, K7));
        MF_MUST(gr_mul(t, ty[3], tw[0], K7));
        MF_MUST(gr_sqr(t, t, K7));
        MF_MUST(gr_sub(rhs, rhs, t, K7));
        _mf_check_equal(lhs, rhs, "addition formula after evictions", tau, K7);

        acb_clear(az); acb_clear(at);
        for (k = 0; k < 4; k++)
            acb_clear(r[k]);
        GR_TMP_CLEAR5(tau, A, P, Q, y, K7);
        GR_TMP_CLEAR4(w, x, t, u, K7);
        GR_TMP_CLEAR4(v, c4, lhs, rhs, K7);
        GR_TMP_CLEAR4(ty[0], ty[1], ty[2], ty[3], K7);
        GR_TMP_CLEAR4(tw[0], tw[1], tw[2], tw[3], K7);
        gr_ctx_clear(K7);
    }

    /* E_k, G_k only for even k >= 2; theta(z) after theta(z/n) at a CM
       point of large degree (the relation not used: independent
       generators, quickly) */
    {
        gr_ctx_t K8;
        gr_ptr tau, z, u, x, th[4];
        acb_t az, at, r[4];
        slong k, n;

        gr_ctx_init_tower_lazy(K8, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_TMP_INIT4(tau, z, u, x, K8);
        GR_TMP_INIT4(th[0], th[1], th[2], th[3], K8);
        acb_init(az); acb_init(at);
        for (k = 0; k < 4; k++)
            acb_init(r[k]);

        MF_MUST(gr_set_str(tau, "1/5 + i", K8));
        for (k = 0; k <= 3; k += 3)
        {
            if (gr_eisenstein_e(x, k, tau, K8) != GR_DOMAIN || gr_eisenstein_g(x, k, tau, K8) != GR_DOMAIN)
            {
                flint_printf("FAIL: E_%wd, G_%wd\n", k, k);
                flint_abort();
            }
        }

        MF_MUST(gr_set_str(tau, "1/5 + 11*i/10", K8));
        MF_MUST(gr_set_str(z, "1/8 + i/4", K8));
        for (n = 3; n <= 6; n += 3)
        {
            MF_MUST(gr_div_ui(u, z, n, K8));
            MF_MUST(gr_jacobi_theta(th[0], th[1], th[2], th[3], u, tau, K8));
            MF_MUST(gr_jacobi_theta(th[0], th[1], th[2], th[3], z, tau, K8));
            MF_MUST(gr_tower_lazy_get_acb(az, z, 128, K8));
            MF_MUST(gr_tower_lazy_get_acb(at, tau, 128, K8));
            acb_modular_theta(r[0], r[1], r[2], r[3], az, at, 128);
            for (k = 0; k < 4; k++)
                _mf_check_overlap(th[k], r[k], "theta_k(z) after theta_k(z/n), large CM", tau, K8);
        }

        acb_clear(az); acb_clear(at);
        for (k = 0; k < 4; k++)
            acb_clear(r[k]);
        GR_TMP_CLEAR4(tau, z, u, x, K8);
        GR_TMP_CLEAR4(th[0], th[1], th[2], th[3], K8);
        gr_ctx_clear(K8);
    }

    gr_ctx_clear(K);
    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
