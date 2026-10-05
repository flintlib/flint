/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "fmpq.h"
#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_elliptic.h"
#include "gr_special.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"

#define SP_PREC 128

/* a random argument: rational, algebraic, complex, or involving pi, log, exp */
static void
_sp_rand_arg(gr_ptr x, flint_rand_t state, int allow_complex, gr_ctx_t K)
{
    gr_ptr t;
    fmpq_t q;
    int kind = n_randint(state, allow_complex ? 7 : 5);

    GR_TMP_INIT(t, K);
    fmpq_init(q);

    fmpq_set_si(q, (slong) n_randint(state, 41) - 20, 1 + n_randint(state, 8));
    GR_MUST_SUCCEED(gr_set_fmpq(x, q, K));

    if (kind == 1)
    {
        /* + sqrt(k) */
        GR_MUST_SUCCEED(gr_set_ui(t, 2 + n_randint(state, 6), K));
        GR_MUST_SUCCEED(gr_sqrt(t, t, K));
        GR_MUST_SUCCEED(gr_add(x, x, t, K));
    }
    else if (kind == 2)
    {
        /* + pi/k */
        GR_MUST_SUCCEED(gr_pi(t, K));
        GR_MUST_SUCCEED(gr_div_ui(t, t, 1 + n_randint(state, 7), K));
        GR_MUST_SUCCEED(gr_add(x, x, t, K));
    }
    else if (kind == 3)
    {
        /* + log(k) */
        GR_MUST_SUCCEED(gr_set_ui(t, 2 + n_randint(state, 5), K));
        GR_MUST_SUCCEED(gr_log(t, t, K));
        GR_MUST_SUCCEED(gr_add(x, x, t, K));
    }
    else if (kind == 4)
    {
        /* + exp(1/k) */
        GR_MUST_SUCCEED(gr_set_ui(t, 1, K));
        GR_MUST_SUCCEED(gr_div_ui(t, t, 1 + n_randint(state, 5), K));
        GR_MUST_SUCCEED(gr_exp(t, t, K));
        GR_MUST_SUCCEED(gr_add(x, x, t, K));
    }
    else if (kind >= 5)
    {
        /* + (r i) */
        GR_MUST_SUCCEED(gr_i(t, K));
        fmpq_set_si(q, (slong) n_randint(state, 21) - 10, 1 + n_randint(state, 4));
        GR_MUST_SUCCEED(gr_mul_fmpq(t, t, q, K));
        GR_MUST_SUCCEED(gr_add(x, x, t, K));
        if (kind == 6)
        {
            GR_MUST_SUCCEED(gr_set_ui(t, 3, K));
            GR_MUST_SUCCEED(gr_sqrt(t, t, K));
            GR_MUST_SUCCEED(gr_add(x, x, t, K));
        }
    }

    fmpq_clear(q);
    GR_TMP_CLEAR(t, K);
}

typedef int (*sp_lazy_func)(gr_ptr, gr_srcptr, gr_ctx_t);

static void
_sp_ref(acb_t res, int which, slong param, const acb_t z, slong prec)
{
    acb_t s;
    fmpz_t k;
    acb_init(s);
    fmpz_init(k);
    switch (which)
    {
        case 0: acb_gamma(res, z, prec); break;
        case 1: acb_hypgeom_erf(res, z, prec); break;
        case 2: acb_hypgeom_erfc(res, z, prec); break;
        case 3: acb_hypgeom_erfi(res, z, prec); break;
        case 4: acb_digamma(res, z, prec); break;
        case 5: acb_set_si(s, param); acb_polygamma(res, s, z, prec); break;
        case 6: fmpz_set_si(k, param); acb_lambertw(res, z, k, 0, prec); break;
        case 7: acb_polylog_si(res, param, z, prec); break;
        case 8: acb_elliptic_k(res, z, prec); break;
        case 9: acb_elliptic_e(res, z, prec); break;
        case 10: acb_rgamma(res, z, prec); break;
        case 11: acb_zeta(res, z, prec); break;
    }
    acb_clear(s);
    fmpz_clear(k);
}

static int
_sp_eval(gr_ptr res, int which, slong param, gr_srcptr x, gr_ctx_t K)
{
    gr_ptr s;
    fmpz_t k;
    int status;
    GR_TMP_INIT(s, K);
    fmpz_init(k);
    GR_MUST_SUCCEED(gr_set_si(s, param, K));
    fmpz_set_si(k, param);
    switch (which)
    {
        case 0: status = gr_gamma(res, x, K); break;
        case 1: status = gr_erf(res, x, K); break;
        case 2: status = gr_erfc(res, x, K); break;
        case 3: status = gr_erfi(res, x, K); break;
        case 4: status = gr_digamma(res, x, K); break;
        case 5: status = gr_polygamma(res, s, x, K); break;
        case 6: status = gr_lambertw_fmpz(res, x, k, K); break;
        case 7: status = gr_polylog(res, s, x, K); break;
        case 8: status = gr_elliptic_k(res, x, K); break;
        case 9: status = gr_elliptic_e(res, x, K); break;
        case 10: status = gr_rgamma(res, x, K); break;
        case 11: status = gr_zeta(res, x, K); break;
        default: status = GR_UNABLE;
    }
    GR_TMP_CLEAR(s, K);
    fmpz_clear(k);
    return status;
}

static void
_sp_check_equal(gr_srcptr a, gr_srcptr b, const char * what, gr_srcptr z, gr_ctx_t K)
{
    truth_t t = gr_equal(a, b, K);
    if (t != T_TRUE)
    {
        flint_printf("FAIL: %s (%d)\n", what, t);
        flint_printf("z = "); gr_println(z, K);
        flint_printf("a = "); gr_println(a, K);
        flint_printf("b = "); gr_println(b, K);
        flint_abort();
    }
}

TEST_FUNCTION_START(gr_tower_special, state)
{
    gr_ctx_t QQ, K, R;
    slong iter;

    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
    gr_tower_lazy_ctx_set_print(K, GR_TOWER_PRINT_NUMERIC | GR_TOWER_PRINT_SYMBOLIC | GR_TOWER_PRINT_DEFS, 10);

    /* special values */
    {
        gr_ptr x, y, z;
        GR_TMP_INIT3(x, y, z, K);

#define SP_EXPECT(expr_str, code, ref_str) \
        if (gr_set_str(x, expr_str, K) != GR_SUCCESS || (code) != GR_SUCCESS || gr_set_str(z, ref_str, K) != GR_SUCCESS) \
        { \
            flint_printf("FAIL: %s at %s (evaluation)\n", #code, expr_str); \
            flint_abort(); \
        } \
        _sp_check_equal(y, z, #code " at " expr_str, x, K);

        SP_EXPECT("5", gr_gamma(y, x, K), "24");
        SP_EXPECT("1/2", gr_gamma(y, x, K), "sqrt(pi)");
        SP_EXPECT("-3/2", gr_gamma(y, x, K), "4*sqrt(pi)/3");
        SP_EXPECT("1/2", gr_rgamma(y, x, K), "1/sqrt(pi)");
        SP_EXPECT("-2", gr_rgamma(y, x, K), "0");
        SP_EXPECT("0", gr_erf(y, x, K), "0");
        SP_EXPECT("0", gr_lambertw(y, x, K), "0");
        SP_EXPECT("-4", gr_zeta(y, x, K), "0");
        SP_EXPECT("-1", gr_zeta(y, x, K), "-1/12");
        SP_EXPECT("0", gr_zeta(y, x, K), "-1/2");
        SP_EXPECT("2", gr_zeta(y, x, K), "pi^2/6");
        SP_EXPECT("6", gr_zeta(y, x, K), "pi^6/945");
        SP_EXPECT("1", gr_digamma(y, x, K), "-euler");
        SP_EXPECT("1/2", gr_digamma(y, x, K), "-euler-2*log(2)");
        SP_EXPECT("1/4", gr_digamma(y, x, K), "-euler-pi/2-3*log(2)");
        SP_EXPECT("1/2", gr_dilog(y, x, K), "pi^2/12-log(2)^2/2");
        SP_EXPECT("1", gr_dilog(y, x, K), "pi^2/6");
        SP_EXPECT("-1", gr_dilog(y, x, K), "-pi^2/12");
        SP_EXPECT("0", gr_elliptic_k(y, x, K), "pi/2");
        SP_EXPECT("1", gr_elliptic_e(y, x, K), "1");

        /* Gamma(1/3) Gamma(2/3) = 2 pi / sqrt(3) */
        GR_MUST_SUCCEED(gr_set_str(x, "1/3", K));
        GR_MUST_SUCCEED(gr_gamma(y, x, K));
        GR_MUST_SUCCEED(gr_set_str(x, "2/3", K));
        GR_MUST_SUCCEED(gr_gamma(z, x, K));
        GR_MUST_SUCCEED(gr_mul(y, y, z, K));
        GR_MUST_SUCCEED(gr_set_str(z, "2*pi/sqrt(3)", K));
        _sp_check_equal(y, z, "Gamma(1/3) Gamma(2/3)", x, K);

        /* psi'(1/4) + psi'(3/4) = 2 pi^2 */
        GR_MUST_SUCCEED(gr_set_ui(z, 1, K));
        GR_MUST_SUCCEED(gr_set_str(x, "1/4", K));
        GR_MUST_SUCCEED(gr_polygamma(y, z, x, K));
        GR_MUST_SUCCEED(gr_set_str(x, "3/4", K));
        GR_MUST_SUCCEED(gr_polygamma(x, z, x, K));
        GR_MUST_SUCCEED(gr_add(y, y, x, K));
        GR_MUST_SUCCEED(gr_set_str(z, "2*pi^2", K));
        _sp_check_equal(y, z, "psi'(1/4) + psi'(3/4)", x, K);

        /* the Hurwitz zeta normal form (phase 5): closed forms */
        {
            gr_ptr m;
            GR_TMP_INIT(m, K);
#define SP_PSI(mm, arg, ref) \
            GR_MUST_SUCCEED(gr_set_ui(m, mm, K)); \
            GR_MUST_SUCCEED(gr_set_str(x, arg, K)); \
            GR_MUST_SUCCEED(gr_polygamma(y, m, x, K)); \
            GR_MUST_SUCCEED(gr_set_str(z, ref, K)); \
            _sp_check_equal(y, z, "polygamma(" #mm ", " arg ") = " ref, x, K);
#define SP_LI(ss, arg, ref) \
            GR_MUST_SUCCEED(gr_set_ui(m, ss, K)); \
            GR_MUST_SUCCEED(gr_set_str(x, arg, K)); \
            GR_MUST_SUCCEED(gr_polylog(y, m, x, K)); \
            GR_MUST_SUCCEED(gr_set_str(z, ref, K)); \
            _sp_check_equal(y, z, "polylog(" #ss ", " arg ") = " ref, x, K);

            SP_PSI(1, "1/4", "pi^2+8*catalan");
            SP_PSI(1, "3/4", "pi^2-8*catalan");
            SP_PSI(1, "5/4", "pi^2+8*catalan-16");
            SP_PSI(2, "1/4", "-2*pi^3-56*zeta(3)");
            SP_PSI(2, "1/3", "-4*pi^3/(3*sqrt(3))-26*zeta(3)");
            SP_PSI(2, "1/6", "-4*sqrt(3)*pi^3-182*zeta(3)");
            SP_PSI(3, "1/2", "pi^4");
            SP_LI(2, "i", "-pi^2/48+i*catalan");
            SP_LI(3, "i", "-3*zeta(3)/32+i*pi^3/32");
            /* the golden ratio values (distribution with n = 2 within
               one orbit; the logarithm comes after the generator) */
            SP_LI(2, "(3-sqrt(5))/2", "pi^2/15-log((sqrt(5)-1)/2)^2");
            SP_LI(2, "(sqrt(5)-1)/2", "pi^2/10-log((sqrt(5)-1)/2)^2");
            SP_LI(2, "(1-sqrt(5))/2", "-pi^2/15+log((sqrt(5)+1)/2)^2/2");
            SP_LI(2, "-(1+sqrt(5))/2", "-pi^2/10-log((sqrt(5)+1)/2)^2");
            /* Li_2(w) + Li_2(w^2) = -pi^2/9, w = exp(2 pi i/3) */
            GR_MUST_SUCCEED(gr_set_ui(m, 2, K));
            GR_MUST_SUCCEED(gr_set_str(x, "exp(2*pi*i/3)", K));
            GR_MUST_SUCCEED(gr_polylog(y, m, x, K));
            GR_MUST_SUCCEED(gr_sqr(x, x, K));
            GR_MUST_SUCCEED(gr_polylog(z, m, x, K));
            GR_MUST_SUCCEED(gr_add(y, y, z, K));
            GR_MUST_SUCCEED(gr_set_str(z, "-pi^2/9", K));
            _sp_check_equal(y, z, "Li_2(w) + Li_2(w^2)", x, K);
            GR_TMP_CLEAR(m, K);
#undef SP_PSI
#undef SP_LI
        }

        /* erfc(100) = 1 - erf(100) = 6.4e-4346: beyond the separation
           precision (no structure theorem decides it) */
        {
            int c;
            GR_MUST_SUCCEED(gr_set_str(x, "erfc(100)", K));
            if (gr_is_zero(x, K) != T_FALSE)
            {
                flint_printf("FAIL: erfc(100) != 0\n");
                flint_abort();
            }
            GR_MUST_SUCCEED(gr_set_str(y, "exp(-10000)/(100*sqrt(pi))", K));
            if (gr_cmp(&c, x, y, K) != GR_SUCCESS || c >= 0)
            {
                flint_printf("FAIL: erfc(100) < exp(-10000)/(100 sqrt(pi))\n");
                flint_abort();
            }
        }

        /* W(-1/e) = -1, and W(1/2) is not rational */
        GR_MUST_SUCCEED(gr_set_si(x, -1, K));
        GR_MUST_SUCCEED(gr_exp(x, x, K));
        GR_MUST_SUCCEED(gr_neg(x, x, K));
        GR_MUST_SUCCEED(gr_lambertw(y, x, K));
        GR_MUST_SUCCEED(gr_set_si(z, -1, K));
        _sp_check_equal(y, z, "W(-1/e)", x, K);

        /* Lambert W in Richardson's algorithm: e^W(z) = z / W(z) */
        SP_EXPECT("2*log(2)", gr_lambertw(y, x, K), "log(2)");
        SP_EXPECT("3*exp(3)", gr_lambertw(y, x, K), "3");
        SP_EXPECT("pi*exp(pi)", gr_lambertw(y, x, K), "pi");
        SP_EXPECT("-log(2)/2", gr_lambertw(y, x, K), "-log(2)");
        {
            fmpz_t k;
            fmpz_init_set_si(k, -1);
            SP_EXPECT("-log(2)/2", gr_lambertw_fmpz(y, x, k, K), "-2*log(2)");
            fmpz_clear(k);
        }
        GR_MUST_SUCCEED(gr_set_str(x, "1", K));
        GR_MUST_SUCCEED(gr_lambertw(y, x, K));
        GR_MUST_SUCCEED(gr_exp(z, y, K));
        GR_MUST_SUCCEED(gr_mul(z, z, y, K));
        GR_MUST_SUCCEED(gr_one(x, K));
        _sp_check_equal(z, x, "W(1) exp(W(1)) = 1", x, K);
        GR_MUST_SUCCEED(gr_log(z, y, K));
        GR_MUST_SUCCEED(gr_add(z, z, y, K));
        if (gr_is_zero(z, K) != T_TRUE)
        {
            flint_printf("FAIL: log(W(1)) + W(1) = 0\n");
            flint_abort();
        }

        /* poles */
        GR_MUST_SUCCEED(gr_set_si(x, -3, K));
        if (gr_gamma(y, x, K) != GR_DOMAIN || gr_digamma(y, x, K) != GR_DOMAIN)
        {
            flint_printf("FAIL: poles\n");
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_one(x, K));
        if (gr_zeta(y, x, K) != GR_DOMAIN || gr_elliptic_k(y, x, K) != GR_DOMAIN)
        {
            flint_printf("FAIL: poles (2)\n");
            flint_abort();
        }

        GR_TMP_CLEAR3(x, y, z, K);
    }

    /* numerical agreement with arb (a GR_UNABLE is allowed for any
       single argument, but most random arguments must be evaluated) */
    {
    slong successes = 0, total = 40 * flint_test_multiplier();
    for (iter = 0; iter < total; iter++)
    {
        gr_ptr x, y;
        acb_t z, v, w;
        int which = n_randint(state, 12), status;
        slong param = 0;

        GR_TMP_INIT2(x, y, K);
        acb_init(z);
        acb_init(v);
        acb_init(w);

        if (which == 5)
            param = 1 + n_randint(state, 3);
        else if (which == 6)
            param = (slong) n_randint(state, 3) - 1;
        else if (which == 7)
            param = (slong) n_randint(state, 7) - 2;

        _sp_rand_arg(x, state, 1, K);
        if ((which == 8 || which == 9 || which == 7) && n_randint(state, 2))
            GR_MUST_SUCCEED(gr_div_ui(x, x, 17, K));

        status = _sp_eval(y, which, param, x, K);
        GR_MUST_SUCCEED(gr_tower_lazy_get_acb(z, x, SP_PREC, K));
        _sp_ref(v, which, param, z, SP_PREC);

        if (status == GR_SUCCESS)
        {
            successes++;
            GR_MUST_SUCCEED(gr_tower_lazy_get_acb(w, y, SP_PREC, K));
            if (!acb_overlaps(v, w))
            {
                flint_printf("FAIL: numerical value (function %d, param %wd)\n", which, param);
                flint_printf("x = "); gr_println(x, K);
                flint_printf("y = "); gr_println(y, K);
                acb_printn(v, 30, 0); flint_printf("\n");
                acb_printn(w, 30, 0); flint_printf("\n");
                flint_abort();
            }
        }
        else if (status == GR_DOMAIN)
        {
            if (acb_is_finite(v) && which != 6 && which != 7)
            {
                flint_printf("FAIL: domain error for a finite value (function %d)\n", which);
                flint_printf("x = "); gr_println(x, K);
                flint_abort();
            }
        }

        GR_TMP_CLEAR2(x, y, K);
        acb_clear(z);
        acb_clear(v);
        acb_clear(w);
    }

    if (total >= 8 && successes < total / 4)
    {
        flint_printf("FAIL: only %wd of %wd random evaluations succeeded\n", successes, total);
        flint_abort();
    }

    }

    /* functional equations hold exactly in the representation */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_ptr z, a, b, c, d;
        int status;

        GR_TMP_INIT5(z, a, b, c, d, K);

        _sp_rand_arg(z, state, 1, K);

        /* Gamma(z + 1) = z Gamma(z) */
        status = gr_gamma(a, z, K);
        if (status == GR_SUCCESS)
        {
            GR_MUST_SUCCEED(gr_add_ui(c, z, 1, K));
            GR_MUST_SUCCEED(gr_gamma(b, c, K));
            GR_MUST_SUCCEED(gr_mul(c, a, z, K));
            _sp_check_equal(b, c, "Gamma(z+1) = z Gamma(z)", z, K);

            /* Gamma(z) Gamma(1 - z) sin(pi z) = pi (z not an integer) */
            GR_MUST_SUCCEED(gr_neg(c, z, K));
            GR_MUST_SUCCEED(gr_add_ui(c, c, 1, K));
            if (gr_gamma(b, c, K) != GR_SUCCESS)
                goto skip_reflection;
            GR_MUST_SUCCEED(gr_mul(b, b, a, K));
            GR_MUST_SUCCEED(gr_pi(c, K));
            GR_MUST_SUCCEED(gr_mul(d, c, z, K));
            GR_MUST_SUCCEED(gr_sin(d, d, K));
            GR_MUST_SUCCEED(gr_mul(b, b, d, K));
            _sp_check_equal(b, c, "reflection formula", z, K);

            /* psi(z + 1) = psi(z) + 1/z */
            GR_MUST_SUCCEED(gr_digamma(a, z, K));
            GR_MUST_SUCCEED(gr_add_ui(c, z, 1, K));
            GR_MUST_SUCCEED(gr_digamma(b, c, K));
            GR_MUST_SUCCEED(gr_inv(c, z, K));
            GR_MUST_SUCCEED(gr_add(c, c, a, K));
            _sp_check_equal(b, c, "psi(z+1) = psi(z) + 1/z", z, K);

            /* psi(1 - z) - psi(z) = pi cot(pi z) */
            GR_MUST_SUCCEED(gr_neg(c, z, K));
            GR_MUST_SUCCEED(gr_add_ui(c, c, 1, K));
            GR_MUST_SUCCEED(gr_digamma(b, c, K));
            GR_MUST_SUCCEED(gr_sub(b, b, a, K));
            GR_MUST_SUCCEED(gr_pi(c, K));
            GR_MUST_SUCCEED(gr_mul(c, c, z, K));
            GR_MUST_SUCCEED(gr_sin(d, c, K));
            GR_MUST_SUCCEED(gr_cos(c, c, K));
            GR_MUST_SUCCEED(gr_div(c, c, d, K));
            GR_MUST_SUCCEED(gr_pi(d, K));
            GR_MUST_SUCCEED(gr_mul(c, c, d, K));
            _sp_check_equal(b, c, "psi(1-z) - psi(z) = pi cot(pi z)", z, K);

skip_reflection:

            /* conj(Gamma(z)) = Gamma(conj(z)) */
            GR_MUST_SUCCEED(gr_gamma(a, z, K));
            GR_MUST_SUCCEED(gr_conj(a, a, K));
            GR_MUST_SUCCEED(gr_conj(c, z, K));
            GR_MUST_SUCCEED(gr_gamma(b, c, K));
            _sp_check_equal(a, b, "conj(Gamma(z)) = Gamma(conj(z))", z, K);
        }

        /* erf(-z) = -erf(z), erfc(z) = 1 - erf(z), erfi(z) = -i erf(i z) */
        if (gr_is_zero(z, K) == T_FALSE)
        {
            GR_MUST_SUCCEED(gr_erf(a, z, K));
            GR_MUST_SUCCEED(gr_neg(c, z, K));
            GR_MUST_SUCCEED(gr_erf(b, c, K));
            GR_MUST_SUCCEED(gr_neg(b, b, K));
            _sp_check_equal(a, b, "erf odd", z, K);
            GR_MUST_SUCCEED(gr_erfc(b, z, K));
            GR_MUST_SUCCEED(gr_add(b, b, a, K));
            GR_MUST_SUCCEED(gr_one(c, K));
            _sp_check_equal(b, c, "erf + erfc", z, K);
            GR_MUST_SUCCEED(gr_i(c, K));
            GR_MUST_SUCCEED(gr_mul(c, c, z, K));
            GR_MUST_SUCCEED(gr_erfi(b, c, K));
            GR_MUST_SUCCEED(gr_i(c, K));
            GR_MUST_SUCCEED(gr_mul(b, b, c, K));
            GR_MUST_SUCCEED(gr_neg(b, b, K));
            _sp_check_equal(a, b, "erf(z) = -i erfi(i z)", z, K);
            GR_MUST_SUCCEED(gr_conj(a, a, K));
            GR_MUST_SUCCEED(gr_conj(c, z, K));
            GR_MUST_SUCCEED(gr_erf(b, c, K));
            _sp_check_equal(a, b, "conj(erf(z)) = erf(conj(z))", z, K);
        }

        GR_TMP_CLEAR5(z, a, b, c, d, K);
    }

    /* Gamma at rationals: Gauss's multiplication theorem holds exactly
       across levels (the normal forms are consistent),
       prod_{j<n} Gamma(a + j/n) = (2 pi)^((n-1)/2) n^(1/2 - n a) Gamma(n a) */
    for (iter = 0; iter < 4 * flint_test_multiplier(); iter++)
    {
        gr_ptr x, y, z, t;
        fmpq_t a, b, e;
        slong n = 2 + n_randint(state, 3), q, j;

        GR_TMP_INIT4(x, y, z, t, K);
        fmpq_init(a);
        fmpq_init(b);
        fmpq_init(e);

        do {
            q = 2 + n_randint(state, 36 / n - 1);
            fmpq_set_si(a, 1 + n_randint(state, q), q);
        } while (fmpz_cmp_ui(fmpq_numref(a), 0) <= 0 || fmpq_cmp_ui(a, 1) >= 0 ||
                 fmpz_cmp_ui(fmpq_denref(a), 36 / n) > 0);

        GR_MUST_SUCCEED(gr_one(y, K));
        for (j = 0; j < n; j++)
        {
            fmpq_set_si(b, j, n);
            fmpq_add(b, b, a);
            GR_MUST_SUCCEED(gr_set_fmpq(x, b, K));
            GR_MUST_SUCCEED(gr_gamma(x, x, K));
            GR_MUST_SUCCEED(gr_mul(y, y, x, K));
        }

        fmpq_mul_si(b, a, n);
        GR_MUST_SUCCEED(gr_set_fmpq(z, b, K));
        GR_MUST_SUCCEED(gr_gamma(z, z, K));
        GR_MUST_SUCCEED(gr_pi(t, K));
        GR_MUST_SUCCEED(gr_mul_ui(t, t, 2, K));
        fmpq_set_si(e, n - 1, 2);
        GR_MUST_SUCCEED(gr_pow_fmpq(t, t, e, K));
        GR_MUST_SUCCEED(gr_mul(z, z, t, K));
        GR_MUST_SUCCEED(gr_set_ui(t, n, K));
        fmpq_set_si(e, 1, 2);
        fmpq_sub(e, e, b);
        GR_MUST_SUCCEED(gr_pow_fmpq(t, t, e, K));
        GR_MUST_SUCCEED(gr_mul(z, z, t, K));

        _sp_check_equal(y, z, "Gauss multiplication", x, K);

        fmpq_clear(a);
        fmpq_clear(b);
        fmpq_clear(e);
        GR_TMP_CLEAR4(x, y, z, t, K);
    }

    /* Gamma(1/6) = 2^(-1/3) (3/pi)^(1/2) Gamma(1/3)^2 */
    {
        gr_ptr x, y;
        GR_TMP_INIT2(x, y, K);
        GR_MUST_SUCCEED(gr_set_str(x, "1/6", K));
        GR_MUST_SUCCEED(gr_gamma(x, x, K));
        GR_MUST_SUCCEED(gr_set_str(y, "gamma(1/3)^2 * sqrt(3/pi) / 2^(1/3)", K));
        _sp_check_equal(x, y, "Gamma(1/6)", x, K);
        GR_TMP_CLEAR2(x, y, K);
    }

    /* beyond the normal form (denominators > 36): the values are
       generators (after reflection), dependent on each other; the zero
       test finds the distribution relations between them (a fresh
       context for each case, to keep the towers small) */
    for (iter = 0; iter < 2 * flint_test_multiplier() + 3; iter++)
    {
        gr_ctx_t K2;
        gr_ptr x, y, z, t;
        fmpq_t a, b, e;
        slong n, q, j;

        gr_ctx_init_tower_lazy(K2, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_TMP_INIT4(x, y, z, t, K2);
        fmpq_init(a);
        fmpq_init(b);
        fmpq_init(e);

        if (iter == 0)
        {
            n = 3;
            fmpq_set_si(a, 13, 40);
        }
        else if (iter <= 2)
        {
            /* levels 111 and 164: roots of 3 (resp. 2) of orders 2 and
               37 (resp. 41) next to roots of unity of orders 3 and 4 */
            q = (iter == 1) ? 37 : 41;
            n = (iter == 1) ? 3 : 4;
            fmpq_set_si(a, 1 + n_randint(state, q - 1), q);
        }
        else
        {
            n = 2;
            q = 37 + n_randint(state, 24);
            fmpq_set_si(a, 1 + n_randint(state, q - 1), q);
            if (fmpz_cmp_ui(fmpq_denref(a), 37) < 0)
                fmpq_set_si(a, 1, 40);
        }

        GR_MUST_SUCCEED(gr_one(y, K2));
        for (j = 0; j < n; j++)
        {
            fmpq_set_si(b, j, n);
            fmpq_add(b, b, a);
            GR_MUST_SUCCEED(gr_set_fmpq(x, b, K2));
            GR_MUST_SUCCEED(gr_gamma(x, x, K2));
            GR_MUST_SUCCEED(gr_mul(y, y, x, K2));
        }

        fmpq_mul_si(b, a, n);
        GR_MUST_SUCCEED(gr_set_fmpq(z, b, K2));
        GR_MUST_SUCCEED(gr_gamma(z, z, K2));
        GR_MUST_SUCCEED(gr_pi(t, K2));
        GR_MUST_SUCCEED(gr_mul_ui(t, t, 2, K2));
        fmpq_set_si(e, n - 1, 2);
        GR_MUST_SUCCEED(gr_pow_fmpq(t, t, e, K2));
        GR_MUST_SUCCEED(gr_mul(z, z, t, K2));
        GR_MUST_SUCCEED(gr_set_ui(t, n, K2));
        fmpq_set_si(e, 1, 2);
        fmpq_sub(e, e, b);
        GR_MUST_SUCCEED(gr_pow_fmpq(t, t, e, K2));
        GR_MUST_SUCCEED(gr_mul(z, z, t, K2));

        _sp_check_equal(y, z, "Gauss multiplication (generators)", x, K2);

        fmpq_clear(a);
        fmpq_clear(b);
        fmpq_clear(e);
        GR_TMP_CLEAR4(x, y, z, t, K2);
        gr_ctx_clear(K2);
    }

    /* Gauss's multiplication theorem at irrational arguments: the values
       lie on a line a w + b, related by the zero test */
    for (iter = 0; iter < 2 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t K2;
        gr_ptr x, y, z, t, w;
        fmpq_t b, e;
        slong n, j;
        const char * ws[] = { "sqrt(2)", "sqrt(3)/2", "pi", "log(2)", "sqrt(2)/3", "(1+sqrt(5))/2", "pi/4+1" };

        gr_ctx_init_tower_lazy(K2, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_TMP_INIT5(x, y, z, t, w, K2);
        fmpq_init(b);
        fmpq_init(e);

        n = 2 + n_randint(state, 2);
        GR_MUST_SUCCEED(gr_set_str(w, ws[n_randint(state, 7)], K2));
        fmpq_set_si(b, (slong) n_randint(state, 5) - 2, 1 + n_randint(state, 3));
        GR_MUST_SUCCEED(gr_add_fmpq(w, w, b, K2));       /* z = w + b */

        /* prod_j Gamma(z + j/n) */
        GR_MUST_SUCCEED(gr_one(y, K2));
        for (j = 0; j < n; j++)
        {
            fmpq_set_si(e, j, n);
            GR_MUST_SUCCEED(gr_add_fmpq(x, w, e, K2));
            GR_MUST_SUCCEED(gr_gamma(x, x, K2));
            GR_MUST_SUCCEED(gr_mul(y, y, x, K2));
        }

        /* (2 pi)^((n-1)/2) n^(1/2 - n z) Gamma(n z) */
        GR_MUST_SUCCEED(gr_mul_ui(z, w, n, K2));
        GR_MUST_SUCCEED(gr_gamma(x, z, K2));
        GR_MUST_SUCCEED(gr_pi(t, K2));
        GR_MUST_SUCCEED(gr_mul_ui(t, t, 2, K2));
        fmpq_set_si(e, n - 1, 2);
        GR_MUST_SUCCEED(gr_pow_fmpq(t, t, e, K2));
        GR_MUST_SUCCEED(gr_mul(x, x, t, K2));
        GR_MUST_SUCCEED(gr_neg(z, z, K2));
        fmpq_set_si(e, 1, 2);
        GR_MUST_SUCCEED(gr_add_fmpq(z, z, e, K2));
        GR_MUST_SUCCEED(gr_set_ui(t, n, K2));
        GR_MUST_SUCCEED(gr_pow(t, t, z, K2));
        GR_MUST_SUCCEED(gr_mul(x, x, t, K2));

        _sp_check_equal(y, x, "Gauss multiplication (irrational arguments)", w, K2);

        /* a wrong constant is detected */
        GR_MUST_SUCCEED(gr_mul_ui(x, x, 2, K2));
        if (gr_equal(y, x, K2) != T_FALSE)
        {
            flint_printf("FAIL: Gauss multiplication (irrational arguments), perturbed\n");
            flint_printf("w = "); gr_println(w, K2);
            flint_abort();
        }

        fmpq_clear(b);
        fmpq_clear(e);
        GR_TMP_CLEAR5(x, y, z, t, w, K2);
        gr_ctx_clear(K2);
    }

    /* the real view: real values at real arguments */
    gr_ctx_init_tower_lazy_view(R, K, GR_TOWER_LAZY_REAL);
    for (iter = 0; iter < 5 * flint_test_multiplier(); iter++)
    {
        gr_ptr x, y;
        int status;

        GR_TMP_INIT2(x, y, R);
        _sp_rand_arg(x, state, 0, R);

        status = gr_gamma(y, x, R);
        if (status != GR_SUCCESS && status != GR_DOMAIN)
        {
            flint_printf("FAIL: real Gamma (%d)\n", status);
            gr_println(x, R);
            flint_abort();
        }
        status = gr_erfi(y, x, R);
        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: real erfi (%d)\n", status);
            gr_println(x, R);
            flint_abort();
        }
        if (gr_tower_lazy_is_real(y, R) != T_TRUE)
        {
            flint_printf("FAIL: erfi not real\n");
            flint_abort();
        }

        GR_TMP_CLEAR2(x, y, R);
    }
    gr_ctx_clear(R);

    /* the multiplication theorem of the polygamma functions,
       sum_{k<n} psi^(m)(x + k/n) = n^(m+1) psi^(m)(n x), at rational x:
       within the normal form (shared context), and beyond it (the zero
       test finds the relations between the generators; fresh contexts),
       with a perturbed counterpart which must be found unequal */
    for (iter = 0; iter < 6 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t K2;
        gr_ptr x, y, z, t, mm;
        fmpq_t c, b;
        slong n, m, p, q, k;
        int beyond = (iter % 3 == 2);
        gr_ctx_struct * KK;

        if (beyond)
        {
            gr_ctx_init_tower_lazy(K2, QQ, GR_TOWER_MERGE_EXPRESS);
            KK = K2;
            n = 2 + n_randint(state, 2);
            /* (levels up to 240: q <= 120 for n = 2, q <= 80 for n = 3) */
            q = 61 + n_randint(state, (n == 2) ? 60 : 20);
        }
        else
        {
            KK = K;
            n = 2 + n_randint(state, 3);
            q = 1 + n_randint(state, 60 / n);
        }
        m = 1 + n_randint(state, 4);
        p = 1 + n_randint(state, 2 * q);

        GR_TMP_INIT5(x, y, z, t, mm, KK);
        fmpq_init(c);
        fmpq_init(b);

        GR_MUST_SUCCEED(gr_set_si(mm, m, KK));
        GR_MUST_SUCCEED(gr_zero(y, KK));
        for (k = 0; k < n; k++)
        {
            fmpq_set_si(c, p, q);
            fmpq_set_si(b, k, n);
            fmpq_add(c, c, b);
            GR_MUST_SUCCEED(gr_set_fmpq(x, c, KK));
            GR_MUST_SUCCEED(gr_polygamma(t, mm, x, KK));
            GR_MUST_SUCCEED(gr_add(y, y, t, KK));
        }
        fmpq_set_si(c, p * n, q);
        GR_MUST_SUCCEED(gr_set_fmpq(x, c, KK));
        GR_MUST_SUCCEED(gr_polygamma(z, mm, x, KK));
        GR_MUST_SUCCEED(gr_set_ui(t, n, KK));
        GR_MUST_SUCCEED(gr_pow_ui(t, t, m + 1, KK));
        GR_MUST_SUCCEED(gr_mul(z, z, t, KK));
        _sp_check_equal(y, z, "polygamma multiplication", x, KK);

        /* perturbed: psi^(m)(x) replaced by psi^(m)(x + 1/(n q)) in the sum */
        if (beyond || iter % 2 == 0)
        {
            fmpq_set_si(c, p, q);
            GR_MUST_SUCCEED(gr_set_fmpq(x, c, KK));
            GR_MUST_SUCCEED(gr_polygamma(t, mm, x, KK));
            GR_MUST_SUCCEED(gr_sub(y, y, t, KK));
            fmpq_set_si(b, 1, n * q);
            fmpq_add(c, c, b);
            GR_MUST_SUCCEED(gr_set_fmpq(x, c, KK));
            GR_MUST_SUCCEED(gr_polygamma(t, mm, x, KK));
            GR_MUST_SUCCEED(gr_add(y, y, t, KK));
            if (gr_equal(y, z, KK) != T_FALSE)
            {
                flint_printf("FAIL: perturbed polygamma multiplication (m = %wd, n = %wd, x = %wd/%wd)\n", m, n, p, q);
                flint_abort();
            }
        }

        fmpq_clear(c);
        fmpq_clear(b);
        GR_TMP_CLEAR5(x, y, z, t, mm, KK);
        if (beyond)
            gr_ctx_clear(K2);
    }

    /* polylogarithms at roots of unity: Li_s(w^2) = 2^(s-1) (Li_s(w) +
       Li_s(-w)), w = exp(2 pi i a/N) */
    for (iter = 0; iter < 3 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t K2;
        gr_ptr w, x, y, z, t, ss;
        /* (sometimes at levels 2N up to 240, in the normal form through
           the cached lattices) */
        slong sv = 2 + n_randint(state, 3), N = (iter % 3 == 2) ? 31 + n_randint(state, 90) : 3 + n_randint(state, 28);
        slong a = 1 + n_randint(state, N - 1);
        char str[64];

        gr_ctx_init_tower_lazy(K2, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_TMP_INIT5(w, x, y, z, t, K2);
        GR_TMP_INIT(ss, K2);

        flint_sprintf(str, "exp(2*pi*i*%wd/%wd)", a, N);
        GR_MUST_SUCCEED(gr_set_str(w, str, K2));
        GR_MUST_SUCCEED(gr_set_si(ss, sv, K2));
        GR_MUST_SUCCEED(gr_sqr(t, w, K2));
        GR_MUST_SUCCEED(gr_polylog(x, ss, t, K2));
        GR_MUST_SUCCEED(gr_polylog(y, ss, w, K2));
        GR_MUST_SUCCEED(gr_neg(t, w, K2));
        GR_MUST_SUCCEED(gr_polylog(z, ss, t, K2));
        GR_MUST_SUCCEED(gr_add(y, y, z, K2));
        GR_MUST_SUCCEED(gr_mul_2exp_si(y, y, sv - 1, K2));
        _sp_check_equal(x, y, "polylog duplication", w, K2);

        GR_TMP_CLEAR5(w, x, y, z, t, K2);
        GR_TMP_CLEAR(ss, K2);
        gr_ctx_clear(K2);
    }

    /* the polygamma multiplication theorem at irrational arguments: the
       values lie on a line a w + b, related by the zero test (shift,
       reflection and multiplication along the line), with a perturbed
       counterpart */
    for (iter = 0; iter < 3 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t K2;
        gr_ptr w, x, y, z, t, mm;
        fmpq_t b, e;
        slong n, m, j;
        int pert;
        const char * ws[] = { "sqrt(2)", "sqrt(3)/2", "pi", "log(2)", "sqrt(2)/3", "(1+sqrt(5))/2", "pi/4+1", "i*sqrt(2)+1/3" };

        gr_ctx_init_tower_lazy(K2, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_TMP_INIT5(w, x, y, z, t, K2);
        GR_TMP_INIT(mm, K2);
        fmpq_init(b);
        fmpq_init(e);

        n = 2 + n_randint(state, 2);
        m = 1 + n_randint(state, 3);
        pert = n_randint(state, 2);
        GR_MUST_SUCCEED(gr_set_str(w, ws[n_randint(state, 8)], K2));
        fmpq_set_si(b, (slong) n_randint(state, 5) - 2, 1 + n_randint(state, 3));
        GR_MUST_SUCCEED(gr_add_fmpq(w, w, b, K2));
        GR_MUST_SUCCEED(gr_set_si(mm, m, K2));

        GR_MUST_SUCCEED(gr_zero(y, K2));
        for (j = 0; j < n; j++)
        {
            fmpq_set_si(e, j, n);
            if (pert && j == 1)
                fmpq_set_si(e, 1, n + 1);
            GR_MUST_SUCCEED(gr_add_fmpq(x, w, e, K2));
            GR_MUST_SUCCEED(gr_polygamma(t, mm, x, K2));
            GR_MUST_SUCCEED(gr_add(y, y, t, K2));
        }
        GR_MUST_SUCCEED(gr_mul_ui(x, w, n, K2));
        GR_MUST_SUCCEED(gr_polygamma(z, mm, x, K2));
        GR_MUST_SUCCEED(gr_set_ui(t, n, K2));
        GR_MUST_SUCCEED(gr_pow_ui(t, t, m + 1, K2));
        GR_MUST_SUCCEED(gr_mul(z, z, t, K2));

        if (!pert)
            _sp_check_equal(y, z, "polygamma multiplication on a line", w, K2);
        else if (gr_equal(y, z, K2) != T_FALSE)
        {
            flint_printf("FAIL: perturbed polygamma multiplication on a line\n");
            gr_println(w, K2);
            flint_abort();
        }

        fmpq_clear(b);
        fmpq_clear(e);
        GR_TMP_CLEAR5(w, x, y, z, t, K2);
        GR_TMP_CLEAR(mm, K2);
        gr_ctx_clear(K2);
    }

    /* the dilogarithm under the anharmonic group: the six values of an
       orbit share one generator, and the reflection and inversion
       formulas hold (z nonreal: no branch cut is met); values agree with
       acb_polylog */
    for (iter = 0; iter < 4 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t K2;
        gr_ptr z, a, b, c, d, two;
        acb_t u, v;
        fmpq_t re, im;
        int which = n_randint(state, 3);

        gr_ctx_init_tower_lazy(K2, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_TMP_INIT5(z, a, b, c, d, K2);
        GR_TMP_INIT(two, K2);
        acb_init(u);
        acb_init(v);
        fmpq_init(re);
        fmpq_init(im);

        fmpq_set_si(re, (slong) n_randint(state, 13) - 6, 1 + n_randint(state, 4));
        fmpq_set_si(im, 1 + n_randint(state, 5), 1 + n_randint(state, 4));
        if (n_randint(state, 2))
            fmpq_neg(im, im);
        GR_MUST_SUCCEED(gr_i(z, K2));
        GR_MUST_SUCCEED(gr_mul_fmpq(z, z, im, K2));
        GR_MUST_SUCCEED(gr_add_fmpq(z, z, re, K2));
        if (n_randint(state, 2))
        {
            /* an irrational real part */
            GR_MUST_SUCCEED(gr_set_ui(a, 2, K2));
            GR_MUST_SUCCEED(gr_sqrt(a, a, K2));
            GR_MUST_SUCCEED(gr_add(z, z, a, K2));
        }
        GR_MUST_SUCCEED(gr_set_ui(two, 2, K2));
        GR_MUST_SUCCEED(gr_dilog(a, z, K2));

        /* the value */
        GR_MUST_SUCCEED(gr_tower_lazy_get_acb(u, z, 128, K2));
        acb_set_si(v, 2);
        acb_polylog(v, v, u, 128);
        GR_MUST_SUCCEED(gr_tower_lazy_get_acb(u, a, 128, K2));
        if (!acb_overlaps(u, v))
        {
            flint_printf("FAIL: dilog value\n");
            gr_println(z, K2);
            gr_println(a, K2);
            flint_abort();
        }

        GR_MUST_SUCCEED(gr_pi(d, K2));
        GR_MUST_SUCCEED(gr_sqr(d, d, K2));
        if (which == 0)
        {
            /* Li_2(z) + Li_2(1 - z) = pi^2/6 - log(z) log(1 - z) */
            GR_MUST_SUCCEED(gr_sub_ui(c, z, 1, K2));
            GR_MUST_SUCCEED(gr_neg(c, c, K2));
            GR_MUST_SUCCEED(gr_dilog(b, c, K2));
            GR_MUST_SUCCEED(gr_add(a, a, b, K2));
            GR_MUST_SUCCEED(gr_log(c, c, K2));
            GR_MUST_SUCCEED(gr_log(b, z, K2));
            GR_MUST_SUCCEED(gr_mul(c, c, b, K2));
            GR_MUST_SUCCEED(gr_div_ui(d, d, 6, K2));
            GR_MUST_SUCCEED(gr_sub(d, d, c, K2));
        }
        else if (which == 1)
        {
            /* Li_2(z) + Li_2(1/z) = -pi^2/6 - log(-z)^2/2 */
            GR_MUST_SUCCEED(gr_inv(c, z, K2));
            GR_MUST_SUCCEED(gr_dilog(b, c, K2));
            GR_MUST_SUCCEED(gr_add(a, a, b, K2));
            GR_MUST_SUCCEED(gr_neg(c, z, K2));
            GR_MUST_SUCCEED(gr_log(c, c, K2));
            GR_MUST_SUCCEED(gr_sqr(c, c, K2));
            GR_MUST_SUCCEED(gr_div_ui(c, c, 2, K2));
            GR_MUST_SUCCEED(gr_div_ui(d, d, 6, K2));
            GR_MUST_SUCCEED(gr_add(d, d, c, K2));
            GR_MUST_SUCCEED(gr_neg(d, d, K2));
        }
        else
        {
            /* Li_2(z) + Li_2(z/(z - 1)) = -log(1 - z)^2/2 */
            GR_MUST_SUCCEED(gr_sub_ui(c, z, 1, K2));
            GR_MUST_SUCCEED(gr_div(c, z, c, K2));
            GR_MUST_SUCCEED(gr_dilog(b, c, K2));
            GR_MUST_SUCCEED(gr_add(a, a, b, K2));
            GR_MUST_SUCCEED(gr_sub_ui(c, z, 1, K2));
            GR_MUST_SUCCEED(gr_neg(c, c, K2));
            GR_MUST_SUCCEED(gr_log(c, c, K2));
            GR_MUST_SUCCEED(gr_sqr(c, c, K2));
            GR_MUST_SUCCEED(gr_div_ui(c, c, 2, K2));
            GR_MUST_SUCCEED(gr_neg(d, c, K2));
        }
        _sp_check_equal(a, d, "dilog functional equation", z, K2);

        fmpq_clear(re);
        fmpq_clear(im);
        acb_clear(u);
        acb_clear(v);
        GR_TMP_CLEAR5(z, a, b, c, d, K2);
        GR_TMP_CLEAR(two, K2);
        gr_ctx_clear(K2);
    }

    /* real arguments: Li_2(1/3) + Li_2(2/3), Li_2(3), Li_2(-1/2) */
    {
        gr_ptr x, y, z;
        GR_TMP_INIT3(x, y, z, K);
        GR_MUST_SUCCEED(gr_set_str(x, "2/3", K));
        GR_MUST_SUCCEED(gr_dilog(y, x, K));
        GR_MUST_SUCCEED(gr_set_str(x, "1/3", K));
        GR_MUST_SUCCEED(gr_dilog(z, x, K));
        GR_MUST_SUCCEED(gr_add(y, y, z, K));
        GR_MUST_SUCCEED(gr_set_str(z, "pi^2/6 - log(1/3)*log(2/3)", K));
        _sp_check_equal(y, z, "Li_2(1/3) + Li_2(2/3)", x, K);
        GR_MUST_SUCCEED(gr_set_str(x, "3", K));
        GR_MUST_SUCCEED(gr_dilog(y, x, K));
        GR_MUST_SUCCEED(gr_set_str(x, "1/3", K));
        GR_MUST_SUCCEED(gr_dilog(z, x, K));
        GR_MUST_SUCCEED(gr_add(y, y, z, K));
        GR_MUST_SUCCEED(gr_set_str(z, "pi^2/3 - log(3)^2/2 - pi*i*log(3)", K));
        _sp_check_equal(y, z, "Li_2(3) + Li_2(1/3)", x, K);
        GR_MUST_SUCCEED(gr_set_str(x, "-1/2", K));
        GR_MUST_SUCCEED(gr_dilog(y, x, K));
        GR_MUST_SUCCEED(gr_set_str(x, "1/3", K));
        GR_MUST_SUCCEED(gr_dilog(z, x, K));
        GR_MUST_SUCCEED(gr_add(y, y, z, K));
        GR_MUST_SUCCEED(gr_set_str(z, "-log(3/2)^2/2", K));
        _sp_check_equal(y, z, "Li_2(-1/2) + Li_2(1/3)", x, K);
        GR_TMP_CLEAR3(x, y, z, K);
    }

    /* complete elliptic integrals: values (imaginary-modulus
       transformation, singular values) against acb; Legendre's relation
       and Landen's transformation found by the zero test, with perturbed
       counterparts */
    for (iter = 0; iter < 4 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t K2;
        gr_ptr m, m1, s, k1, b, ka, ea, kb, eb, t, y;
        acb_t u, v;
        int pert = n_randint(state, 3) == 0;
        char str[100];
        static const char * sv[] = { "3-2*sqrt(2)", "2*sqrt(2)-2", "(2-sqrt(3))/4", "(2+sqrt(3))/4", "17-12*sqrt(2)", "12*sqrt(2)-16", "1/2" };

        gr_ctx_init_tower_lazy(K2, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_TMP_INIT5(m, m1, s, k1, b, K2);
        GR_TMP_INIT5(ka, ea, kb, eb, t, K2);
        GR_TMP_INIT(y, K2);
        acb_init(u);
        acb_init(v);

        if (iter % 4 == 3)
            GR_MUST_SUCCEED(gr_set_str(m, sv[n_randint(state, 7)], K2));
        else if (iter % 2 == 0)
        {
            flint_sprintf(str, "%wd/%wd", 1 + n_randint(state, 8), 9 + n_randint(state, 3));
            GR_MUST_SUCCEED(gr_set_str(m, str, K2));
        }
        else
        {
            flint_sprintf(str, "%wd/%wd + %wd*i/%wd", (slong) n_randint(state, 7) - 3, 1 + n_randint(state, 3),
                (slong) (1 + n_randint(state, 4)) * (n_randint(state, 2) ? 1 : -1), 1 + n_randint(state, 3));
            GR_MUST_SUCCEED(gr_set_str(m, str, K2));
        }

        /* values */
        GR_MUST_SUCCEED(gr_elliptic_k(ka, m, K2));
        GR_MUST_SUCCEED(gr_elliptic_e(ea, m, K2));
        GR_MUST_SUCCEED(gr_tower_lazy_get_acb(u, m, 128, K2));
        acb_elliptic_k(v, u, 128);
        GR_MUST_SUCCEED(gr_tower_lazy_get_acb(u, ka, 128, K2));
        if (!acb_overlaps(u, v))
        {
            flint_printf("FAIL: elliptic K value\n");
            gr_println(m, K2);
            gr_println(ka, K2);
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_tower_lazy_get_acb(u, m, 128, K2));
        acb_elliptic_e(v, u, 128);
        GR_MUST_SUCCEED(gr_tower_lazy_get_acb(u, ea, 128, K2));
        if (!acb_overlaps(u, v))
        {
            flint_printf("FAIL: elliptic E value\n");
            gr_println(m, K2);
            gr_println(ea, K2);
            flint_abort();
        }

        if (iter % 4 != 3)
        {
            /* Legendre: E K' + E' K - K K' = pi/2 */
            GR_MUST_SUCCEED(gr_sub_ui(m1, m, 1, K2));
            GR_MUST_SUCCEED(gr_neg(m1, m1, K2));
            GR_MUST_SUCCEED(gr_elliptic_k(kb, m1, K2));
            GR_MUST_SUCCEED(gr_elliptic_e(eb, m1, K2));
            GR_MUST_SUCCEED(gr_mul(y, ea, kb, K2));
            GR_MUST_SUCCEED(gr_mul(t, eb, ka, K2));
            GR_MUST_SUCCEED(gr_add(y, y, t, K2));
            GR_MUST_SUCCEED(gr_mul(t, ka, kb, K2));
            GR_MUST_SUCCEED(gr_sub(y, y, t, K2));
            GR_MUST_SUCCEED(gr_pi(t, K2));
            GR_MUST_SUCCEED(gr_div_ui(t, t, 2, K2));
            if (pert)
                GR_MUST_SUCCEED(gr_mul_ui(t, t, 2, K2));
            if (gr_equal(y, t, K2) != (pert ? T_FALSE : T_TRUE))
            {
                flint_printf("FAIL: Legendre's relation (perturbed: %d)\n", pert);
                gr_println(m, K2);
                flint_abort();
            }

            /* Landen: K(m) = (1 + k1) K(k1^2), E(m) = (1 + s) E(k1^2) - s K(m) */
            GR_MUST_SUCCEED(gr_sqrt(s, m1, K2));
            GR_MUST_SUCCEED(gr_sub_ui(k1, s, 1, K2));
            GR_MUST_SUCCEED(gr_neg(k1, k1, K2));
            GR_MUST_SUCCEED(gr_add_ui(t, s, 1, K2));
            GR_MUST_SUCCEED(gr_div(k1, k1, t, K2));
            GR_MUST_SUCCEED(gr_sqr(b, k1, K2));
            GR_MUST_SUCCEED(gr_elliptic_k(y, b, K2));
            GR_MUST_SUCCEED(gr_add_ui(t, k1, 1, K2));
            if (pert)
                GR_MUST_SUCCEED(gr_add_ui(t, t, 1, K2));
            GR_MUST_SUCCEED(gr_mul(y, y, t, K2));
            if (gr_equal(ka, y, K2) != (pert ? T_FALSE : T_TRUE))
            {
                flint_printf("FAIL: Landen K (perturbed: %d)\n", pert);
                gr_println(m, K2);
                flint_abort();
            }
            GR_MUST_SUCCEED(gr_elliptic_e(y, b, K2));
            GR_MUST_SUCCEED(gr_add_ui(t, s, 1, K2));
            GR_MUST_SUCCEED(gr_mul(y, y, t, K2));
            GR_MUST_SUCCEED(gr_mul(t, ka, s, K2));
            GR_MUST_SUCCEED(gr_sub(y, y, t, K2));
            if (!pert)
                _sp_check_equal(ea, y, "Landen E", m, K2);
        }

        acb_clear(u);
        acb_clear(v);
        GR_TMP_CLEAR5(m, m1, s, k1, b, K2);
        GR_TMP_CLEAR5(ka, ea, kb, eb, t, K2);
        GR_TMP_CLEAR(y, K2);
        gr_ctx_clear(K2);
    }

    /* the distribution relations of the dilogarithm: Li_2(y^n) =
       n sum_k Li_2(zeta_n^k y) for n = 2, 3, 4 (random rational y, a few
       nonreal and irrational ones), detected by the zero test; perturbed
       versions are nonzero */
    for (iter = 0; iter < 4 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t K2;
        gr_ptr y, x, s, w, t;
        slong n = 2 + n_randint(state, 3), k, kind = n_randint(state, 4);
        int pert = n_randint(state, 4) == 0;
        char str[64];
        fmpq_t q;

        gr_ctx_init_tower_lazy(K2, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_TMP_INIT5(y, x, s, w, t, K2);
        fmpq_init(q);

        /* y: rational (|y| < 1 or > 1), 1/sqrt(2), or (1 + i)/3 */
        do
        {
            fmpq_set_si(q, (slong) n_randint(state, 15) - 7, 1 + n_randint(state, 9));
        }
        while (fmpz_is_zero(fmpq_numref(q)) || fmpz_equal(fmpq_numref(q), fmpq_denref(q)) ||
            (fmpz_is_one(fmpq_denref(q)) && fmpz_equal_si(fmpq_numref(q), -1)));
        if (kind <= 1)
            GR_MUST_SUCCEED(gr_set_fmpq(y, q, K2));
        else if (kind == 2)
            GR_MUST_SUCCEED(gr_set_str(y, n == 2 ? "1/sqrt(2)" : "2/7", K2));
        else
            GR_MUST_SUCCEED(gr_set_str(y, n == 4 ? "(1+i)/3" : "2+i", K2));

        GR_MUST_SUCCEED(gr_pow_ui(x, y, n, K2));
        GR_MUST_SUCCEED(gr_dilog(x, x, K2));
        GR_MUST_SUCCEED(gr_zero(s, K2));
        for (k = 0; k < n; k++)
        {
            flint_sprintf(str, "exp(2*pi*i*%wd/%wd)", k, n);
            GR_MUST_SUCCEED(gr_set_str(w, str, K2));
            GR_MUST_SUCCEED(gr_mul(w, w, y, K2));
            GR_MUST_SUCCEED(gr_dilog(t, w, K2));
            GR_MUST_SUCCEED(gr_add(s, s, t, K2));
        }
        GR_MUST_SUCCEED(gr_mul_ui(s, s, pert ? 2 * n : n, K2));

        if (gr_equal(x, s, K2) != (pert ? T_FALSE : T_TRUE))
        {
            flint_printf("FAIL: dilog distribution (n = %wd, perturbed: %d)\n", n, pert);
            gr_println(y, K2);
            GR_MUST_SUCCEED(gr_sub(x, x, s, K2));
            gr_println(x, K2);
            flint_abort();
        }

        fmpq_clear(q);
        GR_TMP_CLEAR5(y, x, s, w, t, K2);
        gr_ctx_clear(K2);
    }

    /* real forms in a real view: the elementary parts of gamma and
       polygamma values at rationals in the tangent normal form (no roots
       of unity, no i), with the values of the complex field */
    {
        static const slong cases[][3] = { {-1, 3, 5}, {-1, 5, 7}, {2, 3, 5}, {0, 1, 7}, {1, 16, 17}, {1, 4, 5}, {0, 2, 9} };
        gr_ctx_t C, R;
        gr_ptr x, y, z, s;
        slong j;
        char * str, * p;

        gr_ctx_init_tower_lazy(C, QQ, GR_TOWER_MERGE_EXPRESS);
        gr_ctx_init_tower_lazy_view(R, C, GR_TOWER_LAZY_REAL);
        x = gr_heap_init(R); y = gr_heap_init(R); z = gr_heap_init(C); s = gr_heap_init(R);

        for (j = 0; j < (slong) (sizeof(cases) / sizeof(cases[0])); j++)
        {
            fmpq_t c;
            fmpq_init(c);
            fmpq_set_si(c, cases[j][1], cases[j][2]);
            GR_MUST_SUCCEED(gr_set_fmpq(x, c, R));
            GR_MUST_SUCCEED(gr_set_si(s, cases[j][0], R));
            if (cases[j][0] < 0)
                GR_MUST_SUCCEED(gr_gamma(y, x, R));
            else
                GR_MUST_SUCCEED(gr_polygamma(y, s, x, R));
            GR_MUST_SUCCEED(gr_get_str(&str, y, R));
            for (p = str; *p; p++)
                if (p[0] == 'p' && p[1] == 'i')
                    p[0] = p[1] = 'P';
            if (strstr(str, "exp(") != NULL || strstr(str, "i") != NULL)
            {
                flint_printf("FAIL: real form (%wd, %wd/%wd): %s\n", cases[j][0], cases[j][1], cases[j][2], str);
                flint_abort();
            }
            flint_free(str);

            GR_MUST_SUCCEED(gr_set_fmpq(z, c, C));
            if (cases[j][0] < 0)
                GR_MUST_SUCCEED(gr_gamma(z, z, C));
            else
            {
                gr_ptr sc = gr_heap_init(C);
                GR_MUST_SUCCEED(gr_set_si(sc, cases[j][0], C));
                GR_MUST_SUCCEED(gr_polygamma(z, sc, z, C));
                gr_heap_clear(sc, C);
            }
            if (gr_equal(y, z, C) != T_TRUE)
            {
                flint_printf("FAIL: real form vs complex value (%wd, %wd/%wd)\n", cases[j][0], cases[j][1], cases[j][2]);
                flint_abort();
            }
            fmpq_clear(c);
        }

        gr_heap_clear(x, R); gr_heap_clear(y, R); gr_heap_clear(z, C); gr_heap_clear(s, R);
        gr_ctx_clear(R);
        gr_ctx_clear(C);
    }

    /* integer real parts with a nonzero imaginary part */
    {
        gr_ptr x, y, mm;
        GR_TMP_INIT3(x, y, mm, K);
        GR_MUST_SUCCEED(gr_set_str(x, "2+i*sqrt(2)", K));
        GR_MUST_SUCCEED(gr_set_si(mm, 3, K));
        GR_MUST_SUCCEED(gr_polygamma(y, mm, x, K));
        GR_MUST_SUCCEED(gr_gamma(y, x, K));
        GR_TMP_CLEAR3(x, y, mm, K);
    }

    gr_ctx_clear(K);
    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
