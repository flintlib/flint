/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Generated using Claude Code (Fable 5, Opus 5.5).

   Not checked in full detail. */

/* Tests of acb_ode_fundamental_matrix against reference values, mostly
   adapted from doctests of the ore_algebra package for SageMath. */

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_mat.h"
#include "acb_ode.h"
#include "acb_poly.h"
#include "gr.h"
#include "gr_ore_poly.h"
#include "gr_poly.h"

/****************************************************************************/

/* Return zero unless got overlaps ref and has a relative accuracy of at
   least min_acc bits. When ref contains zero, we instead require got to be
   bounded by 2^-min_acc in absolute value. */
static int
check_entry(const acb_t got, const acb_t ref, slong min_acc)
{
    mag_t m;
    int ok;

    FLINT_ASSERT(acb_is_finite(ref));

    if (!acb_overlaps(got, ref))
        return 0;
    if (!acb_contains_zero(ref))
        return acb_rel_accuracy_bits(got) >= min_acc;

    FLINT_ASSERT(mag_cmp_2exp_si(arb_radref(acb_realref(ref)), -30) <= 0);
    FLINT_ASSERT(mag_cmp_2exp_si(arb_radref(acb_imagref(ref)), -30) <= 0);

    mag_init(m);
    acb_get_mag(m, got);
    ok = mag_cmp_2exp_si(m, -min_acc) <= 0;
    mag_clear(m);

    return ok;
}

static void
check_acb(const acb_t got, const acb_t ref, slong min_acc, const char * name)
{
    if (!check_entry(got, ref, min_acc))
    {
        flint_printf("FAIL: %s\n", name);
        flint_printf("got = %{acb}\n", got);
        flint_printf("ref = %{acb}\n", ref);
        flint_printf("accuracy = %wd bits, required %wd\n",
                     acb_rel_accuracy_bits(got), min_acc);
        flint_abort();
    }
}

static void
check_status(int status, const char * name)
{
    if (status != GR_SUCCESS)
    {
        flint_printf("FAIL (unexpected failure): %s\n", name);
        flint_abort();
    }
}

static void
set_acb_str(acb_t res, const char * s, slong prec)
{
    gr_ctx_t CC;

    gr_ctx_init_complex_acb(CC, prec);
    GR_MUST_SUCCEED(gr_set_str(res, s, CC));
    gr_ctx_clear(CC);
}

/****************************************************************************
 * Ordinary points
 ****************************************************************************/

static int
my_exp(acb_t res, const acb_t z, slong prec)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_mat_t mat;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_mat_init(mat, 1, 1);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "D - 1", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, z,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    acb_swap(res, acb_mat_entry(mat, 0, 0));

    acb_mat_clear(mat);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);

    return status;
}

/* sage: (Dx - 1).numerical_solution(ini=[1], path=[0, 1], eps=1e-50)
   [2.7182818284590452353602874713526624977572470936999...]
   [ore_algebra] */
static void
test_exp(void)
{
    acb_t z, res, ref;
    slong prec = 192;

    acb_init(z);
    acb_init(res);
    acb_init(ref);

    acb_one(z);
    check_status(my_exp(res, z, prec), "exp");
    acb_exp(ref, z, prec);
    check_acb(res, ref, prec - 16, "exp");

    acb_clear(ref);
    acb_clear(res);
    acb_clear(z);
}

/* sage: (Dx - 1).numerical_solution([1], [0, i + pi])
   [12.5029695888765...] + [19.4722214188416...]*I
   [ore_algebra] */
static void
test_exp_complex(void)
{
    acb_t z, res, ref;
    slong prec = 128;

    acb_init(z);
    acb_init(res);
    acb_init(ref);

    set_acb_str(z, "pi + I", prec);
    check_status(my_exp(res, z, prec), "exp_complex");
    acb_exp(ref, z, prec);
    check_acb(res, ref, prec - 16, "exp_complex");

    acb_clear(ref);
    acb_clear(res);
    acb_clear(z);
}

/* sage: stdfun.cos.dop.numerical_solution(stdfun.cos.ini,
   ....:                                   [stdfun.cos.pt, pi/3])
   [0.500000000000000...]
   sage: stdfun.sin.dop.numerical_solution(stdfun.sin.ini,
   ....:                                   [stdfun.sin.pt, pi/3])
   [0.86602540378443...]
   with stdfun.cos.dop = stdfun.sin.dop = Dx^2 + 1.
   [ore_algebra] */
static void
test_cos_sin(void)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_t pt;
    acb_poly_t x, s, c;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_init(pt);
    acb_poly_init(x);
    acb_poly_init(s);
    acb_poly_init(c);
    acb_mat_init(mat, 2, 2);

    set_acb_str(pt, "pi/3", prec);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "D^2 + 1", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "cos_sin");

    acb_poly_set_coeff_si(x, 1, 1);
    acb_poly_set_coeff_acb(x, 0, pt);
    acb_poly_sin_cos_series(s, c, x, 2, prec);
    for (slong i = 0; i < 2; i++)
    {
        check_acb(acb_mat_entry(mat, i, 0), c->coeffs + i, prec - 16,
                  "cos_sin");
        check_acb(acb_mat_entry(mat, i, 1), s->coeffs + i, prec - 16,
                  "cos_sin");
    }

    acb_mat_clear(mat);
    acb_poly_clear(c);
    acb_poly_clear(s);
    acb_poly_clear(x);
    acb_clear(pt);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* sage: ((x + 1)*Dx^2 + Dx).numerical_transition_matrix([0,1/2])
   [ [1.00000000000000...] [0.405465108108164...]]
   [             [+/- ...] [0.666666666666666...]]
   [ore_algebra] */
static void
test_log1p(void)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_t pt, c;
    acb_poly_t x, f;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_init(pt);
    acb_init(c);
    acb_poly_init(x);
    acb_poly_init(f);
    acb_mat_init(mat, 2, 2);

    acb_set_d(pt, 0.5);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "(x + 1)*D^2 + D", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "log1p");

    acb_one(c);
    check_acb(acb_mat_entry(mat, 0, 0), c, prec - 16, "log1p");
    acb_zero(c);
    check_acb(acb_mat_entry(mat, 1, 0), c, prec - 16, "log1p");

    acb_poly_set_coeff_si(x, 1, 1);
    acb_poly_set_coeff_acb(x, 0, pt);
    acb_poly_log1p_series(f, x, 2, prec);
    for (slong i = 0; i < 2; i++)
        check_acb(acb_mat_entry(mat, i, 1), f->coeffs + i, prec - 16,
                  "log1p");

    acb_mat_clear(mat);
    acb_poly_clear(f);
    acb_poly_clear(x);
    acb_clear(c);
    acb_clear(pt);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* sage: ((x + 1)*Dx^3 + Dx^2).numerical_transition_matrix([0,1/2])
   [ [1.000000000000000...] [0.5000000000000000...] [0.2163953243244931...]]
   [              [+/- ...]  [1.000000000000000...] [0.8109302162163287...]]
   [              [+/- ...]               [+/- ...] [0.6666666666666666...]]
   The third basis element is 2*((1 + x)*log(1 + x) - x).
   [ore_algebra] */
static void
test_log1p_order3(void)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_t pt, c;
    acb_poly_t x, f, g;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_init(pt);
    acb_init(c);
    acb_poly_init(x);
    acb_poly_init(f);
    acb_poly_init(g);
    acb_mat_init(mat, 3, 3);

    acb_set_d(pt, 0.5);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "(x + 1)*D^3 + D^2", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "log1p_order3");

    /* 1, x */
    for (slong i = 0; i < 3; i++)
    {
        for (slong j = 0; j < 2; j++)
        {
            if (i == 0 && j == 1)
                acb_set(c, pt);
            else
                acb_set_si(c, i == j);
            check_acb(acb_mat_entry(mat, i, j), c, prec - 16, "log1p_order3");
        }
    }

    acb_poly_set_coeff_si(x, 1, 1);
    acb_poly_set_coeff_acb(x, 0, pt);
    acb_poly_log1p_series(f, x, 3, prec);
    acb_poly_add_si(g, x, 1, prec);
    acb_poly_mullow(f, f, g, 3, prec);
    acb_poly_sub(f, f, x, prec);
    acb_poly_scalar_mul_2exp_si(f, f, 1);
    for (slong i = 0; i < 3; i++)
        check_acb(acb_mat_entry(mat, i, 2), f->coeffs + i, prec - 16,
                  "log1p_order3");

    acb_mat_clear(mat);
    acb_poly_clear(g);
    acb_poly_clear(f);
    acb_poly_clear(x);
    acb_clear(c);
    acb_clear(pt);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* sage: ((x + 1)*Dx^3 + Dx).numerical_transition_matrix([0,1/2])
   [  [1.000000000000000...]  [0.4815453970799961...]  [0.2456596136789682...]]
   [               [+/- ...]  [0.8936357901691244...]  [0.9667328760004665...]]
   [               [+/- ...] [-0.1959698689702905...]  [0.9070244207738327...]]
   [ore_algebra] */
static void
test_ordinary_order3(void)
{
    const char * ref_str[] = {
        "1", "0.4815453970799961 +/- 1e-16", "0.2456596136789682 +/- 1e-16",
        "0", "0.8936357901691244 +/- 1e-16", "0.9667328760004665 +/- 1e-16",
        "0", "-0.1959698689702905 +/- 1e-16", "0.9070244207738327 +/- 1e-16",
    };
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_t pt, ref;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_init(pt);
    acb_init(ref);
    acb_mat_init(mat, 3, 3);

    acb_set_d(pt, 0.5);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "(x + 1)*D^3 + D", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "ordinary_order3");

    for (slong i = 0; i < 3; i++)
    {
        for (slong j = 0; j < 3; j++)
        {
            set_acb_str(ref, ref_str[3 * i + j], prec);
            check_acb(acb_mat_entry(mat, i, j), ref, prec - 16,
                      "ordinary_order3");
        }
    }

    acb_mat_clear(mat);
    acb_clear(ref);
    acb_clear(pt);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* Dawson function, as the element of leading monomial x of the local basis
   at 0 of Dx^2 + 2*x*Dx + 2 */
static int
my_dawson(acb_t res, const acb_t z, slong prec)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_mat_t mat;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_mat_init(mat, 1, 2);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "D^2 + 2*x*D + 2", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, z,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    acb_swap(res, acb_mat_entry(mat, 0, 1));

    acb_mat_clear(mat);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);

    return status;
}

/* Dawson function (sqrt(pi)/2)*exp(-z^2)*erfi(z) */
static void
ref_dawson(acb_t res, const acb_t z, slong prec)
{
    acb_t c;

    acb_init(c);

    acb_hypgeom_erfi(res, z, prec);
    acb_sqr(c, z, prec);
    acb_neg(c, c);
    acb_exp(c, c, prec);
    acb_mul(res, res, c, prec);
    set_acb_str(c, "sqrt(pi)/2", prec);
    acb_mul(res, res, c, prec);

    acb_clear(c);
}

/* sage: stdfun.dawson.dop.numerical_solution(stdfun.dawson.ini,
   ....:                                      [stdfun.dawson.pt, 1])
   [0.538079506912768...]
   with stdfun.dawson.dop = Dx^2 + 2*x*Dx + 2.
   [ore_algebra] */
static void
test_dawson(void)
{
    acb_t z, res, ref;
    slong prec = 128;

    acb_init(z);
    acb_init(res);
    acb_init(ref);

    acb_one(z);
    check_status(my_dawson(res, z, prec), "dawson");
    ref_dawson(ref, z, prec);
    check_acb(res, ref, prec - 16, "dawson");

    acb_clear(ref);
    acb_clear(res);
    acb_clear(z);
}

/* erf(z) as a combination of the local basis at 0 of Dx^2 + 2*x*Dx */
static int
my_erf(acb_t res, const acb_t z, slong prec)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_mat_t mat;
    acb_t c;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_mat_init(mat, 1, 2);
    acb_init(c);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "D^2 + 2*x*D", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, z,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    set_acb_str(c, "2/sqrt(pi)", prec);
    acb_mul(res, acb_mat_entry(mat, 0, 1), c, prec);

    acb_clear(c);
    acb_mat_clear(mat);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);

    return status;
}

/* sage: stdfun.erf.dop.numerical_solution(stdfun.erf.ini,
   ....:                                   [stdfun.erf.pt, 1])
   [0.842700792949714...]
   with stdfun.erf.dop = Dx^2 + 2*x*Dx, stdfun.erf.ini = [0, 2/sqrt(pi)].
   [ore_algebra] */
static void
test_erf(void)
{
    acb_t z, res, ref;
    slong prec = 128;

    acb_init(z);
    acb_init(res);
    acb_init(ref);

    acb_one(z);
    check_status(my_erf(res, z, prec), "erf");
    acb_hypgeom_erf(ref, z, prec);
    check_acb(res, ref, prec - 16, "erf");

    acb_clear(ref);
    acb_clear(res);
    acb_clear(z);
}

/* sage: (Dx^2 + 2*x*Dx).numerical_solution([0, 2/sqrt(pi)], [0,i])
   [+/- ...] + [1.65042575879754...]*I
   [ore_algebra] */
static void
test_erf_imag(void)
{
    acb_t z, res, ref;
    slong prec = 128;

    acb_init(z);
    acb_init(res);
    acb_init(ref);

    acb_onei(z);
    check_status(my_erf(res, z, prec), "erf_imag");
    acb_hypgeom_erf(ref, z, prec);
    check_acb(res, ref, prec - 16, "erf_imag");

    acb_clear(ref);
    acb_clear(res);
    acb_clear(z);
}

/* erfi(z) as a combination of the local basis at 0 of Dx^2 - 2*x*Dx */
static int
my_erfi(acb_t res, const acb_t z, slong prec)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_mat_t mat;
    acb_t c;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_mat_init(mat, 1, 2);
    acb_init(c);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "D^2 - 2*x*D", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, z,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    set_acb_str(c, "2/sqrt(pi)", prec);
    acb_mul(res, acb_mat_entry(mat, 0, 1), c, prec);

    acb_clear(c);
    acb_mat_clear(mat);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);

    return status;
}

/* sage: stdfun.erfi.dop.numerical_solution(stdfun.erfi.ini,
   ....:                                    [stdfun.erfi.pt, 1])
   [1.65042575879754...]
   with stdfun.erfi.dop = Dx^2 - 2*x*Dx, stdfun.erfi.ini = [0, 2/sqrt(pi)].
   [ore_algebra] */
static void
test_erfi(void)
{
    acb_t z, res, ref;
    slong prec = 128;

    acb_init(z);
    acb_init(res);
    acb_init(ref);

    acb_one(z);
    check_status(my_erfi(res, z, prec), "erfi");
    acb_hypgeom_erfi(ref, z, prec);
    check_acb(res, ref, prec - 16, "erfi");

    acb_clear(ref);
    acb_clear(res);
    acb_clear(z);
}

/* Set res to the jet of length len of Ai at z */
static int
my_airy_ai_jet(acb_ptr res, const acb_t z, slong len, slong prec)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_mat_t mat;
    acb_t zero, ai0, aip0;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_mat_init(mat, len, 2);
    acb_init(zero);
    acb_init(ai0);
    acb_init(aip0);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "D^2 - x", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, z,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    acb_hypgeom_airy(ai0, aip0, NULL, NULL, zero, prec);
    for (slong i = 0; i < len; i++)
    {
        acb_mul(res + i, acb_mat_entry(mat, i, 0), ai0, prec);
        acb_addmul(res + i, acb_mat_entry(mat, i, 1), aip0, prec);
    }

    acb_clear(aip0);
    acb_clear(ai0);
    acb_clear(zero);
    acb_mat_clear(mat);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);

    return status;
}

/* sage: stdfun.airy_ai.dop.numerical_solution(stdfun.airy_ai.ini,
   ....:                                       [stdfun.airy_ai.pt, 1])
   [0.135292416312881...]
   with stdfun.airy_ai.dop = Dx^2 - x.
   [ore_algebra] */
static void
test_airy_ai(void)
{
    acb_t z, res, ref;
    slong prec = 128;

    acb_init(z);
    acb_init(res);
    acb_init(ref);

    acb_one(z);
    check_status(my_airy_ai_jet(res, z, 1, prec), "airy_ai");
    acb_hypgeom_airy(ref, NULL, NULL, NULL, z, prec);
    check_acb(res, ref, prec - 16, "airy_ai");

    acb_clear(ref);
    acb_clear(res);
    acb_clear(z);
}

/* sage: (Dx^2-x).numerical_solution(ini, [0,2], post_transform=Dx^10)
   [2.34553207877...]
   sage: airy_ai(10, 2.)
   2.345532078777...
   with ini = [Ai(0), Ai'(0)].
   [ore_algebra] */
static void
test_airy_ai_jet(void)
{
    acb_t z;
    acb_ptr res, ref;
    slong prec = 128;

    acb_init(z);
    res = _acb_vec_init(11);
    ref = _acb_vec_init(11);

    acb_set_si(z, 2);
    check_status(my_airy_ai_jet(res, z, 11, prec), "airy_ai_jet");
    acb_hypgeom_airy_jet(ref, NULL, z, 11, prec);
    for (slong i = 0; i < 11; i++)
        check_acb(res + i, ref + i, prec - 32, "airy_ai_jet");

    _acb_vec_clear(ref, 11);
    _acb_vec_clear(res, 11);
    acb_clear(z);
}

/* sage: (Dx - i).numerical_transition_matrix([0, sqrt(2)])
   [[0.1559436947653...] + [0.9877659459927...]*I]
   Operator over ZZ[i].
   [ore_algebra] */
static void
test_exp_gaussian(void)
{
    gr_ctx_t ZZi, Pol, Dop;
    gr_ore_poly_t dop;
    acb_t pt, ref;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpzi(ZZi);
    gr_ctx_init_gr_poly(Pol, ZZi);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_init(pt);
    acb_init(ref);
    acb_mat_init(mat, 1, 1);

    set_acb_str(pt, "sqrt(2)", prec);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "D - I", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "exp_gaussian");

    acb_mul_onei(ref, pt);
    acb_exp(ref, ref, prec);
    check_acb(acb_mat_entry(mat, 0, 0), ref, prec - 16, "exp_gaussian");

    acb_mat_clear(mat);
    acb_clear(ref);
    acb_clear(pt);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZi);
}

/* sage: (Dx - zeta7).numerical_transition_matrix([0, 1])
   [[1.32375209616333...] + [1.31434281345999...]*I]
   with zeta7 = exp(2*pi*i/7). Operator over QQbar.
   [ore_algebra] */
static void
test_exp_zeta7(void)
{
    gr_ctx_t QQbar, Pol, Dop;
    gr_ore_poly_t dop;
    gr_ptr zeta7;
    acb_t pt, ref;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_complex_qqbar(QQbar);
    gr_ctx_init_gr_poly(Pol, QQbar);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    zeta7 = gr_heap_init(QQbar);
    acb_init(pt);
    acb_init(ref);
    acb_mat_init(mat, 1, 1);

    /* gr_ore_poly_set_str does not support "D - (-1)^(2/7)" */
    acb_one(pt);
    GR_MUST_SUCCEED(gr_set_str(zeta7, "(-1)^(2/7)", QQbar));
    GR_MUST_SUCCEED(gr_ore_poly_gen(dop, Dop));
    GR_MUST_SUCCEED(gr_ore_poly_sub_other(dop, dop, zeta7, QQbar, Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "exp_zeta7");

    set_acb_str(ref, "exp(exp(2*pi*I/7))", prec);
    check_acb(acb_mat_entry(mat, 0, 0), ref, prec - 16, "exp_zeta7");

    acb_mat_clear(mat);
    acb_clear(ref);
    acb_clear(pt);
    gr_heap_clear(zeta7, QQbar);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(QQbar);
}

/* sage: (Du^2 - RBF(pi)).numerical_solution(ini=[1,0], path=[0,1])
   [3.02759638974277...]
   Operator with ball coefficients; the basis is cosh(s*x), sinh(s*x)/s
   where s = sqrt(pi).
   [ore_algebra] */
static void
test_cosh_sqrt_pi(void)
{
    gr_ctx_t CC, Pol, Dop;
    gr_ore_poly_t dop;
    acb_t pt, s;
    acb_poly_t x, sh, ch;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_complex_acb(CC, prec);
    gr_ctx_init_gr_poly(Pol, CC);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_init(pt);
    acb_init(s);
    acb_poly_init(x);
    acb_poly_init(sh);
    acb_poly_init(ch);
    acb_mat_init(mat, 2, 2);

    acb_one(pt);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "D^2 - pi", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "cosh_sqrt_pi");

    set_acb_str(s, "sqrt(pi)", prec);
    acb_poly_set_coeff_si(x, 1, 1);
    acb_poly_set_coeff_acb(x, 0, pt);
    acb_poly_scalar_mul(x, x, s, prec);
    acb_poly_sinh_cosh_series(sh, ch, x, 2, prec);
    acb_poly_scalar_div(sh, sh, s, prec);
    for (slong i = 0; i < 2; i++)
    {
        check_acb(acb_mat_entry(mat, i, 0), ch->coeffs + i, prec - 16,
                  "cosh_sqrt_pi");
        check_acb(acb_mat_entry(mat, i, 1), sh->coeffs + i, prec - 16,
                  "cosh_sqrt_pi");
    }

    acb_mat_clear(mat);
    acb_poly_clear(ch);
    acb_poly_clear(sh);
    acb_poly_clear(x);
    acb_clear(s);
    acb_clear(pt);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(CC);
}

/* sage: (Du - CBF(i*pi)).numerical_transition_matrix([0,1])
   [[-1.0000000000000...] + [+/- ...]*I]
   [ore_algebra] */
static void
test_exp_i_pi(void)
{
    gr_ctx_t CC, Pol, Dop;
    gr_ore_poly_t dop;
    acb_t pt, ref;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_complex_acb(CC, prec);
    gr_ctx_init_gr_poly(Pol, CC);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_init(pt);
    acb_init(ref);
    acb_mat_init(mat, 1, 1);

    acb_one(pt);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "D - pi*I", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "exp_i_pi");

    acb_set_si(ref, -1);
    check_acb(acb_mat_entry(mat, 0, 0), ref, prec - 16, "exp_i_pi");

    acb_mat_clear(mat);
    acb_clear(ref);
    acb_clear(pt);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(CC);
}

/* Solution of Neher's IVP (see below) */
static int
neher(acb_t res, const acb_t z, slong prec)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_mat_t mat;
    acb_t c;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_mat_init(mat, 1, 4);
    acb_init(c);

    /* the operator is scaled by 2 to clear denominators */
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop,
            "2*D^4 - (2*x^2 + 20*x + 52)*D^3 + (40*x + 199)*D^2"
            " - (2*x^2 + 20*x + 50)*D + 4*x^2 + 8*x - 59", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, z,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);

    acb_mul_si(res, acb_mat_entry(mat, 0, 0), 5, prec);
    acb_addmul_si(res, acb_mat_entry(mat, 0, 1), 4, prec);
    acb_set_d(c, 1.5);
    acb_addmul(res, acb_mat_entry(mat, 0, 2), c, prec);
    acb_div_si(c, acb_mat_entry(mat, 0, 3), 3, prec);
    acb_add(res, res, c, prec);

    acb_clear(c);
    acb_mat_clear(mat);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);

    return status;
}

/* sage: dop = (Dx^4 - (x^2+10*x+26)*Dx^3 - (-20*x-99-1/2)*Dx^2
   ....:        - (x^2+10*x+25)*Dx - (-2*x^2-4*x+29+1/2))
   sage: ini = [5,4,3/2,1/3]
   sage: dop.numerical_solution(ini, [0,1])
   [10.873127314 +/- ...e-10]
   (Example from M. Neher, 1999; known to be subject to precision loss.)
   [ore_algebra] */
static void
test_neher(void)
{
    acb_t z, res, ref;
    slong prec = 128;

    acb_init(z);
    acb_init(res);
    acb_init(ref);

    acb_one(z);
    check_status(neher(res, z, prec), "neher");
    set_acb_str(ref, "10.873127314 +/- 1e-9", prec);
    check_acb(res, ref, prec - 64, "neher");

    acb_clear(ref);
    acb_clear(res);
    acb_clear(z);
}

/* sage: dop.numerical_solution(ini, [0,3/2], 1e-30, two_point_mode=False)
   [15.68591174618322...]
   (continued from test_neher)
   [ore_algebra] */
static void
test_neher_3_2(void)
{
    acb_t z, res, ref;
    slong prec = 128;

    acb_init(z);
    acb_init(res);
    acb_init(ref);

    acb_set_d(z, 1.5);
    check_status(neher(res, z, prec), "neher_3_2");
    set_acb_str(ref, "15.68591174618322 +/- 1e-14", prec);
    check_acb(res, ref, prec - 64, "neher_3_2");

    acb_clear(ref);
    acb_clear(res);
    acb_clear(z);
}

/* sage: (((x^200-1)//(x-1))*Dx^2 + 5*x*Dx
   ....:   + sum((i+3)*x**i for i in range(500))).numerical_transition_matrix(
   ....:                                                    [0,1/16])
   [ [0.99412284058692...] [0.062180605822178...]]
   [[-0.18799810365157...]  [0.98478282225787...]]
   [ore_algebra] */
static void
test_high_degree(void)
{
    const char * ref_str[] = {
        "0.99412284058692 +/- 1e-14", "0.062180605822178 +/- 1e-15",
        "-0.18799810365157 +/- 1e-14", "0.98478282225787 +/- 1e-14",
    };
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_t pt, ref;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_init(pt);
    acb_init(ref);
    acb_mat_init(mat, 2, 2);

    acb_set_d(pt, 0.0625);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop,
            "((x^200 - 1)/(x - 1))*D^2 + 5*x*D"
            " + 3*(x^500 - 1)/(x - 1)"
            " + x*(499*x^500 - 500*x^499 + 1)/(x - 1)^2", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "high_degree");

    for (slong i = 0; i < 2; i++)
    {
        for (slong j = 0; j < 2; j++)
        {
            set_acb_str(ref, ref_str[2 * i + j], prec);
            check_acb(acb_mat_entry(mat, i, j), ref, prec - 16,
                      "high_degree");
        }
    }

    acb_mat_clear(mat);
    acb_clear(ref);
    acb_clear(pt);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* sage: dop = ((x-1)^2*Dx-1).lclm(Dx-1)
   sage: dop.local_basis_monomials(0)
   [1, x^2]
   sage: dop.numerical_solution([1, 1/2],[0, 2], assume_analytic=True)
   [7.389056098930650...] + [+/- ...]*I
   The origin is an apparent singularity with exponents 0, 2. We evaluate
   at 1/2 instead of continuing to 2 through the singular point 1. The lclm,
   which kills exp(x) and exp(x/(1-x)), was computed by hand; the local basis
   at 0 is (3*exp(x) - exp(x/(1-x)))/2, exp(x/(1-x)) - exp(x).
   [ore_algebra] */
static void
test_apparent_singularity(void)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_t pt;
    acb_poly_t x, e, f, g;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_init(pt);
    acb_poly_init(x);
    acb_poly_init(e);
    acb_poly_init(f);
    acb_poly_init(g);
    acb_mat_init(mat, 2, 2);

    acb_set_d(pt, 0.5);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop,
            "x*(x - 2)*(x - 1)^2*D^2 - ((x - 1)^4 + 2*x - 3)*D + x^2 - 2",
            Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "apparent_singularity");

    acb_poly_set_coeff_si(x, 1, 1);
    acb_poly_set_coeff_acb(x, 0, pt);
    acb_poly_exp_series(e, x, 2, prec);
    acb_poly_neg(g, x);
    acb_poly_add_si(g, g, 1, prec);
    acb_poly_div_series(g, x, g, 2, prec);
    acb_poly_exp_series(g, g, 2, prec);
    acb_poly_scalar_mul_2exp_si(f, e, 1);
    acb_poly_add(f, f, e, prec);
    acb_poly_sub(f, f, g, prec);
    acb_poly_scalar_mul_2exp_si(f, f, -1);
    acb_poly_sub(g, g, e, prec);
    for (slong i = 0; i < 2; i++)
    {
        check_acb(acb_mat_entry(mat, i, 0), f->coeffs + i, prec - 16,
                  "apparent_singularity");
        check_acb(acb_mat_entry(mat, i, 1), g->coeffs + i, prec - 16,
                  "apparent_singularity");
    }

    acb_mat_clear(mat);
    acb_poly_clear(g);
    acb_poly_clear(f);
    acb_poly_clear(e);
    acb_poly_clear(x);
    acb_clear(pt);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/****************************************************************************
 * Regular singular points
 ****************************************************************************/

/* sage: (x*Dx + 1).numerical_transition_matrix([0,2], algorithm="binsplit")
   [0.50000000000000000]
   The recurrence has order zero.
   [ore_algebra] */
static void
test_inverse(void)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_t pt;
    acb_poly_t x, f;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_init(pt);
    acb_poly_init(x);
    acb_poly_init(f);
    acb_mat_init(mat, 2, 1);

    acb_set_si(pt, 2);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "x*D + 1", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "inverse");

    acb_poly_set_coeff_si(x, 1, 1);
    acb_poly_set_coeff_acb(x, 0, pt);
    acb_poly_inv_series(f, x, 2, prec);
    for (slong i = 0; i < 2; i++)
        check_acb(acb_mat_entry(mat, i, 0), f->coeffs + i, prec - 16,
                  "inverse");

    acb_mat_clear(mat);
    acb_poly_clear(f);
    acb_poly_clear(x);
    acb_clear(pt);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* sage: (x^2*Dx^2 + x*Dx + 1).numerical_transition_matrix([0, 1/2])
   [ [0.769238901363972...] + [0.638961276313634...]*I
     [0.769238901363972...] + [-0.6389612763136...]*I]
   [...]
   Exponents -i, i; the basis is x^(-i), x^i.
   [ore_algebra] */
static void
test_euler_pm_i(void)
{
    gr_ctx_t ZZ, Pol, Dop, Pol1, Euler, CC;
    gr_ore_poly_t dop, ldop;
    acb_ode_exponents_t expos;
    acb_t pt, expo;
    acb_poly_t x, f;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;
    ulong flags = 0;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    acb_ode_local_diffop_ctx_init(Pol1, Euler, Dop, ZZ);
    gr_ctx_init_complex_acb(CC, prec);
    gr_ore_poly_init(dop, Dop);
    gr_ore_poly_init(ldop, Euler);
    acb_ode_exponents_init(expos);
    acb_init(pt);
    acb_init(expo);
    acb_poly_init(x);
    acb_poly_init(f);
    acb_mat_init(mat, 2, 2);

    acb_set_d(pt, 0.5);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "x^2*D^2 + x*D + 1", Dop));
    status |= acb_ode_local_diffop(ldop, &flags, dop, Dop, Euler);
    status |= acb_ode_exponents(expos, ldop, Euler, CC);
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, expos, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "euler_pm_i");

    acb_poly_set_coeff_si(x, 1, 1);
    acb_poly_set_coeff_acb(x, 0, pt);
    for (slong j = 0; j < 2; j++)
    {
        acb_onei(expo);
        if (arb_is_negative(acb_imagref(expos->groups[j].leader)))
            acb_neg(expo, expo);
        acb_poly_pow_acb_series(f, x, expo, 2, prec);
        for (slong i = 0; i < 2; i++)
            check_acb(acb_mat_entry(mat, i, j), f->coeffs + i, prec - 16,
                      "euler_pm_i");
    }

    acb_mat_clear(mat);
    acb_poly_clear(f);
    acb_poly_clear(x);
    acb_clear(expo);
    acb_clear(pt);
    acb_ode_exponents_clear(expos);
    gr_ore_poly_clear(ldop, Euler);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(CC);
    gr_ctx_clear(Euler);
    gr_ctx_clear(Pol1);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* sage: dop = (x*Dx-1001/2).symmetric_product(Dx-1)
   sage: dop = dop._normalize_base_ring()[-1]
   sage: ref = exp(CBF(1/2))/RBF(2)^(1001/2)
   sage: ref.overlaps(dop.numerical_transition_matrix([0, 1/2], 1e-10)[0,0])
   True
   The normalized operator is 2*x*Dx - 2*x - 1001; the solution is
   x^(1001/2)*exp(x).
   [ore_algebra] */
static void
test_large_exponent(void)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_t pt, ref;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_init(pt);
    acb_init(ref);
    acb_mat_init(mat, 1, 1);

    acb_set_d(pt, 0.5);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "2*x*D - 2*x - 1001", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "large_exponent");

    set_acb_str(ref, "exp(1/2)/2^(1001/2)", prec);
    check_acb(acb_mat_entry(mat, 0, 0), ref, prec - 32, "large_exponent");

    acb_mat_clear(mat);
    acb_clear(ref);
    acb_clear(pt);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* sage: dop = (x*Dx+1001/2).symmetric_product(Dx-1)
   sage: dop = dop._normalize_base_ring()[-1]
   sage: val = dop.numerical_transition_matrix([0, 1/2], 1e-10)[0,0]
   sage: ref = CBF(1/2)^(-1001/2)*exp(CBF(1/2))
   sage: ref.overlaps(val)
   True
   The normalized operator is 2*x*Dx - 2*x + 1001.
   [ore_algebra] */
static void
test_large_negative_exponent(void)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_t pt, ref;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_init(pt);
    acb_init(ref);
    acb_mat_init(mat, 1, 1);

    acb_set_d(pt, 0.5);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "2*x*D - 2*x + 1001", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "large_negative_exponent");

    set_acb_str(ref, "exp(1/2)*2^(1001/2)", prec);
    check_acb(acb_mat_entry(mat, 0, 0), ref, prec - 32,
              "large_negative_exponent");

    acb_mat_clear(mat);
    acb_clear(ref);
    acb_clear(pt);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* The next two operators have order 3 and exponents 0, 1000, 1000 at the
   origin. Their local bases there are e(x)*(1 - c*x^1000),
   e(x)*x^1000*log(x), e(x)*x^1000, where e(x) = exp(x) and c = 1/1000! for
   the first one, and e(x) = 1/(1 + x) and c = 1 for the second one. */

static int
large_dispersion(acb_mat_t mat, const char * dop_str, const acb_t pt,
                 slong prec)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, dop_str, Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);

    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);

    return status;
}

/* Jets of length 3 at pt of the three basis elements */
static void
ref_large_dispersion(acb_poly_struct * ref, const acb_poly_t e, const acb_t c,
                     const acb_t pt, slong prec)
{
    acb_poly_t x;

    acb_poly_init(x);

    acb_poly_set_coeff_si(x, 1, 1);
    acb_poly_set_coeff_acb(x, 0, pt);
    acb_poly_pow_ui_trunc_binexp(ref + 2, x, 1000, 3, prec);
    acb_poly_mullow(ref + 2, ref + 2, e, 3, prec);
    acb_poly_log_series(ref + 1, x, 3, prec);
    acb_poly_mullow(ref + 1, ref + 1, ref + 2, 3, prec);
    acb_poly_scalar_mul(ref + 0, ref + 2, c, prec);
    acb_poly_sub(ref + 0, e, ref + 0, prec);

    acb_poly_clear(x);
}

/* sage: h = CBF(1/2)
   sage: #dop = (Dx-1).lclm(x^2*Dx^2 - x*(2*x+1999)*Dx + (x^2 + 1999*x + 1000^2))
   sage: dop = x^2*Dx^3 + (-3*x^2 - 1997*x)*Dx^2 + (3*x^2 + 3994*x + 998001)*Dx - x^2 - 1997*x - 998001
   sage: mat = dop.numerical_transition_matrix([0,1/2], 1e-5)
   sage: mat[0,0].overlaps(exp(h))
   True
   sage: mat[0,1].overlaps(exp(h)*h^1000*log(h))
   True
   sage: mat[0,2].overlaps(exp(h)*h^1000)
   True

   The tail bounds for n < 1000 are astronomical, as the staircase bound
   includes the values of the rational function near the resonance n = 1000,
   where |ind(n)/n^3| ~ (n/1000 - 1)^2 is tiny. They only become usable well
   past n = 1000. In addition, the coefficients of the first basis element
   lose all accuracy as n approaches 1000, and spurious log terms appear
   unless the working precision is large enough (prec >= ~1536). Hence the
   first column is only checked at high precision.
   [ore_algebra] */
static void
test_large_dispersion_exp(void)
{
    const char * dop_str = "x^2*D^3 + (-3*x^2 - 1997*x)*D^2"
                           " + (3*x^2 + 3994*x + 998001)*D"
                           " - x^2 - 1997*x - 998001";
    const slong precs[] = { 128, 2048 };
    const slong min_acc[] = { 112, 2000 };
    acb_t pt, c;
    acb_poly_t x, e;
    acb_poly_struct * ref;
    acb_mat_t mat;

    acb_init(pt);
    acb_init(c);
    acb_poly_init(x);
    acb_poly_init(e);
    ref = _acb_poly_vec_init(3);
    acb_mat_init(mat, 3, 3);

    acb_set_d(pt, 0.5);

    for (slong k = 0; k < 2; k++)
    {
        slong prec = precs[k];

        check_status(large_dispersion(mat, dop_str, pt, prec),
                     "large_dispersion_exp");

        acb_poly_set_coeff_si(x, 1, 1);
        acb_poly_set_coeff_acb(x, 0, pt);
        acb_poly_exp_series(e, x, 3, prec);
        acb_set_si(c, 1001);
        acb_rgamma(c, c, prec);
        ref_large_dispersion(ref, e, c, pt, prec);

        for (slong i = 0; i < 3; i++)
            for (slong j = (k == 0); j < 3; j++)
                check_acb(acb_mat_entry(mat, i, j), ref[j].coeffs + i,
                          min_acc[k], "large_dispersion_exp");
    }

    acb_mat_clear(mat);
    _acb_poly_vec_clear(ref, 3);
    acb_poly_clear(e);
    acb_poly_clear(x);
    acb_clear(c);
    acb_clear(pt);
}

/* sage: dop = (x^3 + x^2)*Dx^3 + (-1994*x^2 - 1997*x)*Dx^2 + (994007*x + 998001)*Dx + 998001
   sage: mat = dop.numerical_transition_matrix([0, 1/2], 1e-5)
   sage: mat[0,0].overlaps(1/(1+h))
   True
   sage: mat[0,1].overlaps(h^1000/(1+h)*log(h))
   True
   sage: mat[0,2].overlaps(h^1000/(1+h))
   True
   [ore_algebra] */
static void
test_large_dispersion_inv(void)
{
    acb_t pt, c;
    acb_poly_t x, e;
    acb_poly_struct * ref;
    acb_mat_t mat;
    slong prec = 128;

    acb_init(pt);
    acb_init(c);
    acb_poly_init(x);
    acb_poly_init(e);
    ref = _acb_poly_vec_init(3);
    acb_mat_init(mat, 3, 3);

    acb_set_d(pt, 0.5);
    check_status(large_dispersion(mat,
                     "(x^3 + x^2)*D^3 + (-1994*x^2 - 1997*x)*D^2"
                     " + (994007*x + 998001)*D + 998001", pt, prec),
                 "large_dispersion_inv");

    acb_poly_set_coeff_si(x, 1, 1);
    acb_add_si(c, pt, 1, prec);
    acb_poly_set_coeff_acb(x, 0, c);
    acb_poly_inv_series(e, x, 3, prec);
    acb_one(c);
    ref_large_dispersion(ref, e, c, pt, prec);

    for (slong i = 0; i < 3; i++)
        for (slong j = 0; j < 3; j++)
            check_acb(acb_mat_entry(mat, i, j), ref[j].coeffs + i, prec - 16,
                      "large_dispersion_inv");

    acb_mat_clear(mat);
    _acb_poly_vec_clear(ref, 3);
    acb_poly_clear(e);
    acb_poly_clear(x);
    acb_clear(c);
    acb_clear(pt);
}

/* sage: alg = QQbar(-20)^(1/3)
   sage: (x*Dx^2 + Dx + x).numerical_transition_matrix([0, alg], 1e-8)
   [ [3.7849872...] +  [1.7263190...]*I  [1.3140884...] + [-2.3112610...]*I]
   [ [1.0831414...] + [-3.3595150...]*I  [-2.0854436...] + [-0.7923237...]*I]
   Bessel operator; the basis is f1 = log(x) + O(x^2), f2 = J_0(x), so that
   f1 = (pi/2)*Y_0(x) - (euler_gamma - log(2))*J_0(x).
   [ore_algebra] */
static void
test_bessel(void)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_t pt, nu, j0, j1, y0, y1, c, d, ref;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_init(pt);
    acb_init(nu);
    acb_init(j0);
    acb_init(j1);
    acb_init(y0);
    acb_init(y1);
    acb_init(c);
    acb_init(d);
    acb_init(ref);
    acb_mat_init(mat, 2, 2);

    set_acb_str(pt, "(-20)^(1/3)", prec);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "x*D^2 + D + x", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "bessel");

    acb_hypgeom_bessel_j(j0, nu, pt, prec);
    acb_hypgeom_bessel_y(y0, nu, pt, prec);
    acb_one(nu);
    acb_hypgeom_bessel_j(j1, nu, pt, prec);
    acb_hypgeom_bessel_y(y1, nu, pt, prec);
    arb_const_euler(acb_realref(c), prec);
    acb_set_si(d, 2);
    acb_log(d, d, prec);
    acb_sub(c, c, d, prec);
    set_acb_str(d, "pi/2", prec);

    /* f1 = d*Y0 - c*J0, f1' = -d*Y1 + c*J1, f2 = J0, f2' = -J1 */
    acb_mul(ref, d, y0, prec);
    acb_submul(ref, c, j0, prec);
    check_acb(acb_mat_entry(mat, 0, 0), ref, prec - 32, "bessel");
    acb_mul(ref, c, j1, prec);
    acb_submul(ref, d, y1, prec);
    check_acb(acb_mat_entry(mat, 1, 0), ref, prec - 32, "bessel");
    check_acb(acb_mat_entry(mat, 0, 1), j0, prec - 32, "bessel");
    acb_neg(ref, j1);
    check_acb(acb_mat_entry(mat, 1, 1), ref, prec - 32, "bessel");

    acb_mat_clear(mat);
    acb_clear(ref);
    acb_clear(d);
    acb_clear(c);
    acb_clear(y1);
    acb_clear(y0);
    acb_clear(j1);
    acb_clear(j0);
    acb_clear(nu);
    acb_clear(pt);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* sage: ((x/1)*Dx^2 - 1).local_basis_monomials(0)
   [1, x]
   sage: ((x/1)*Dx^2 - 1).numerical_transition_matrix([0, 1])
   [[0.0340875989376363...]   [1.59063685463732...]]
   [[-0.579827135138349...]   [2.27958530233606...]]
   With w = 2*sqrt(x), the basis is
   f1 = w*K_1(w) - (2*euler_gamma - 1)*sqrt(x)*I_1(w) = 1 + x*log(x) + O(x^2)
   and f2 = sqrt(x)*I_1(w), with f1' = -2*K_0(w) - (2*euler_gamma - 1)*I_0(w)
   and f2' = I_0(w).
   [ore_algebra] */
static void
test_bessel_i_k(void)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_t pt, w, nu, i0, i1, k0, k1, c, ref;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_init(pt);
    acb_init(w);
    acb_init(nu);
    acb_init(i0);
    acb_init(i1);
    acb_init(k0);
    acb_init(k1);
    acb_init(c);
    acb_init(ref);
    acb_mat_init(mat, 2, 2);

    acb_one(pt);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "x*D^2 - 1", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "bessel_i_k");

    acb_set_si(w, 2);
    acb_hypgeom_bessel_i(i0, nu, w, prec);
    acb_hypgeom_bessel_k(k0, nu, w, prec);
    acb_one(nu);
    acb_hypgeom_bessel_i(i1, nu, w, prec);
    acb_hypgeom_bessel_k(k1, nu, w, prec);
    arb_const_euler(acb_realref(c), prec);
    acb_mul_2exp_si(c, c, 1);
    acb_sub_ui(c, c, 1, prec);

    acb_mul_2exp_si(ref, k1, 1);
    acb_submul(ref, c, i1, prec);
    check_acb(acb_mat_entry(mat, 0, 0), ref, prec - 32, "bessel_i_k");
    acb_mul_2exp_si(ref, k0, 1);
    acb_addmul(ref, c, i0, prec);
    acb_neg(ref, ref);
    check_acb(acb_mat_entry(mat, 1, 0), ref, prec - 32, "bessel_i_k");
    check_acb(acb_mat_entry(mat, 0, 1), i1, prec - 32, "bessel_i_k");
    check_acb(acb_mat_entry(mat, 1, 1), i0, prec - 32, "bessel_i_k");

    acb_mat_clear(mat);
    acb_clear(ref);
    acb_clear(c);
    acb_clear(k1);
    acb_clear(k0);
    acb_clear(i1);
    acb_clear(i0);
    acb_clear(nu);
    acb_clear(w);
    acb_clear(pt);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* sage: (x*Dx^2 + x + 1).numerical_transition_matrix([0, 1], 1e-10)
   [[-0.006152006884...]    [0.4653635461...]]
   [   [-2.148107776...]  [-0.05672090833...]]
   [ore_algebra] */
static void
test_singular_hardcoded(void)
{
    const char * ref_str[] = {
        "-0.006152006884 +/- 1e-12", "0.4653635461 +/- 1e-10",
        "-2.148107776 +/- 1e-9", "-0.05672090833 +/- 1e-11",
    };
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_t pt, ref;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_init(pt);
    acb_init(ref);
    acb_mat_init(mat, 2, 2);

    acb_one(pt);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "x*D^2 + x + 1", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "singular_hardcoded");

    for (slong i = 0; i < 2; i++)
    {
        for (slong j = 0; j < 2; j++)
        {
            set_acb_str(ref, ref_str[2 * i + j], prec);
            check_acb(acb_mat_entry(mat, i, j), ref, prec - 32,
                      "singular_hardcoded");
        }
    }

    acb_mat_clear(mat);
    acb_clear(ref);
    acb_clear(pt);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* Cosine integral, as the solution of x*Dx^3 + 2*Dx^2 + x*Dx with
   ini = [1, euler_gamma, 0] in the local basis with leading monomials
   log(x), 1, x */
static int
my_ci(acb_t res, const acb_t z, slong prec)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_mat_t mat;
    acb_t c;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_mat_init(mat, 1, 3);
    acb_init(c);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "x*D^3 + 2*D^2 + x*D", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, z,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    arb_const_euler(acb_realref(c), prec);
    acb_mul(res, acb_mat_entry(mat, 0, 1), c, prec);
    acb_add(res, res, acb_mat_entry(mat, 0, 0), prec);

    acb_clear(c);
    acb_mat_clear(mat);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);

    return status;
}

/* sage: dop = x*Dx^3 + 2*Dx^2 + x*Dx
   sage: ini = [1, CBF(euler_gamma), 0]
   sage: dop.numerical_solution(ini, path=[0, sqrt(2)])
   [0.46365280236686...]
   [ore_algebra] */
static void
test_ci(void)
{
    acb_t z, res, ref;
    slong prec = 128;

    acb_init(z);
    acb_init(res);
    acb_init(ref);

    set_acb_str(z, "sqrt(2)", prec);
    check_status(my_ci(res, z, prec), "ci");
    acb_hypgeom_ci(ref, z, prec);
    check_acb(res, ref, prec - 32, "ci");

    acb_clear(ref);
    acb_clear(res);
    acb_clear(z);
}

/* sage: dop.numerical_solution(ini, path=[0, 456/123*i+1])
   [6.1267878728616...] + [-3.39197789100074...]*I
   (continued from test_ci)
   [ore_algebra] */
static void
test_ci_complex(void)
{
    acb_t z, res, ref;
    slong prec = 128;

    acb_init(z);
    acb_init(res);
    acb_init(ref);

    set_acb_str(z, "1 + 456/123*I", prec);
    check_status(my_ci(res, z, prec), "ci_complex");
    acb_hypgeom_ci(ref, z, prec);
    check_acb(res, ref, prec - 32, "ci_complex");

    acb_clear(ref);
    acb_clear(res);
    acb_clear(z);
}

/* sage: dop = (Dx*(x*Dx)^2).lclm(Dx-1)
   sage: dop.numerical_solution([2, 0, 0, 0], [0, 1/4])
   [1.921812055672805...]
   Here dop is replaced by the equivalent operator
   ((x^2 + 3*x + 1)*T^4 - (x^3 + 6*x^2 + 7*x + 1)*T^3)/x^2, T = x*Dx,
   written in terms of Dx. The local exponents at 0 are 0, 0, 0, 1 and the
   basis is log(x)^2/2, log(x), 1, exp(x) - 1.
   [ore_algebra] */
static void
test_log2_lclm(void)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_t pt, c;
    acb_poly_t x, f, g;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_init(pt);
    acb_init(c);
    acb_poly_init(x);
    acb_poly_init(f);
    acb_poly_init(g);
    acb_mat_init(mat, 4, 4);

    acb_set_d(pt, 0.25);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop,
            "x^2*(x^2 + 3*x + 1)*D^4 - x*(x^3 - 11*x - 5)*D^3"
            " - (3*x^3 + 11*x^2 - 4)*D^2 - (x + 1)*(x + 4)*D", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "log2_lclm");

    acb_poly_set_coeff_si(x, 1, 1);
    acb_poly_set_coeff_acb(x, 0, pt);
    acb_poly_log_series(f, x, 4, prec);
    acb_poly_mullow(g, f, f, 4, prec);
    acb_poly_scalar_mul_2exp_si(g, g, -1);
    for (slong i = 0; i < 4; i++)
    {
        check_acb(acb_mat_entry(mat, i, 0), g->coeffs + i, prec - 32,
                  "log2_lclm");
        check_acb(acb_mat_entry(mat, i, 1), f->coeffs + i, prec - 32,
                  "log2_lclm");
        acb_set_si(c, i == 0);
        check_acb(acb_mat_entry(mat, i, 2), c, prec - 32, "log2_lclm");
    }
    acb_poly_exp_series(f, x, 4, prec);
    acb_poly_add_si(f, f, -1, prec);
    for (slong i = 0; i < 4; i++)
        check_acb(acb_mat_entry(mat, i, 3), f->coeffs + i, prec - 32,
                  "log2_lclm");

    acb_mat_clear(mat);
    acb_poly_clear(g);
    acb_poly_clear(f);
    acb_poly_clear(x);
    acb_clear(c);
    acb_clear(pt);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* Set res to the jet of length len at z of x^expo * pFq(a; b; x) */
static void
xpow_pfq_jet(acb_ptr res, const acb_t expo,
             acb_srcptr a, slong p, acb_srcptr b, slong q,
             const acb_t z, slong len, slong prec)
{
    acb_poly_struct * ap = _acb_poly_vec_init(p);
    acb_poly_struct * bp = _acb_poly_vec_init(q + 1);
    acb_poly_t x, f, g;

    acb_poly_init(x);
    acb_poly_init(f);
    acb_poly_init(g);

    for (slong i = 0; i < p; i++)
        acb_poly_set_acb(ap + i, a + i);
    for (slong i = 0; i < q; i++)
        acb_poly_set_acb(bp + i, b + i);
    /* the series computed by acb_hypgeom_pfq_series_direct lacks the 1/k! */
    acb_poly_one(bp + q);

    acb_poly_set_coeff_si(x, 1, 1);
    acb_poly_set_coeff_acb(x, 0, z);
    acb_hypgeom_pfq_series_direct(f, ap, p, bp, q + 1, x, 0, -1, len, prec);
    acb_poly_pow_acb_series(g, x, expo, len, prec);
    acb_poly_mullow(g, g, f, len, prec);
    for (slong i = 0; i < len; i++)
        acb_poly_get_coeff_acb(res + i, g, i);

    acb_poly_clear(g);
    acb_poly_clear(f);
    acb_poly_clear(x);
    _acb_poly_vec_clear(bp, q + 1);
    _acb_poly_vec_clear(ap, p);
}

/* Set res to the jet of length len at z of the Whittaker function
   M_{kappa,mu}(x) = exp(-x/2)*x^(mu+1/2)*1F1(mu-kappa+1/2; 1+2*mu; x) */
static void
whittaker_m_jet(acb_ptr res, const acb_t kappa, const acb_t mu,
                const acb_t z, slong len, slong prec)
{
    acb_t expo, a, b;
    acb_ptr h, f, g;

    acb_init(expo);
    acb_init(a);
    acb_init(b);
    h = _acb_vec_init(2);
    f = _acb_vec_init(len);
    g = _acb_vec_init(len);

    acb_set_d(expo, 0.5);
    acb_add(expo, expo, mu, prec);
    acb_sub(a, expo, kappa, prec);
    acb_mul_2exp_si(b, mu, 1);
    acb_add_si(b, b, 1, prec);
    xpow_pfq_jet(f, expo, a, 1, b, 1, z, len, prec);

    /* g = jet of exp(-x/2) */
    acb_mul_2exp_si(h, z, -1);
    acb_neg(h, h);
    acb_set_d(h + 1, -0.5);
    _acb_poly_exp_series(g, h, 2, len, prec);

    _acb_poly_mullow(res, f, len, g, len, len, prec);

    _acb_vec_clear(g, len);
    _acb_vec_clear(f, len);
    _acb_vec_clear(h, 2);
    acb_clear(b);
    acb_clear(a);
    acb_clear(expo);
}

/* sage: dop = 4*x^2*Dx^2 + (-x^2+8*x-11)
   sage: dop.numerical_transition_matrix([0, 10])
   [[-3.829367993175840...]  [7.857756823216673...]]
   [[-1.135875563239369...]  [1.426170676718429...]]
   Exponents 1/2 -/+ sqrt(3); the basis consists of the Whittaker functions
   M_{2,-sqrt(3)} and M_{2,sqrt(3)}.
   [ore_algebra] */
static void
test_whittaker(void)
{
    gr_ctx_t ZZ, Pol, Dop, Pol1, Euler, CC;
    gr_ore_poly_t dop, ldop;
    acb_ode_exponents_t expos;
    acb_t pt, kappa, mu;
    acb_ptr f;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;
    ulong flags = 0;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    acb_ode_local_diffop_ctx_init(Pol1, Euler, Dop, ZZ);
    gr_ctx_init_complex_acb(CC, prec);
    gr_ore_poly_init(dop, Dop);
    gr_ore_poly_init(ldop, Euler);
    acb_ode_exponents_init(expos);
    acb_init(pt);
    acb_init(kappa);
    acb_init(mu);
    f = _acb_vec_init(2);
    acb_mat_init(mat, 2, 2);

    acb_set_si(pt, 10);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "4*x^2*D^2 - x^2 + 8*x - 11",
                                        Dop));
    status |= acb_ode_local_diffop(ldop, &flags, dop, Dop, Euler);
    status |= acb_ode_exponents(expos, ldop, Euler, CC);
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, expos, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "whittaker");

    acb_set_si(kappa, 2);
    for (slong j = 0; j < 2; j++)
    {
        acb_set_si(mu, 3);
        acb_sqrt(mu, mu, prec);
        if (arb_is_negative(acb_realref(expos->groups[j].leader)))
            acb_neg(mu, mu);
        whittaker_m_jet(f, kappa, mu, pt, 2, prec);
        for (slong i = 0; i < 2; i++)
            check_acb(acb_mat_entry(mat, i, j), f + i, prec - 32,
                      "whittaker");
    }

    acb_mat_clear(mat);
    _acb_vec_clear(f, 2);
    acb_clear(mu);
    acb_clear(kappa);
    acb_clear(pt);
    acb_ode_exponents_clear(expos);
    gr_ore_poly_clear(ldop, Euler);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(CC);
    gr_ctx_clear(Euler);
    gr_ctx_clear(Pol1);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* sage: ((x*Dx - 3/7)).lclm(Dx - 1).numerical_transition_matrix([0,1/3])
   [[1.395612425086089...] + [+/- ...]*I  [0.6244813348581596...]]
   [...]
   The lclm was computed by hand. The basis is exp(x), x^(3/7).
   [ore_algebra] */
static void
test_lclm_xpow(void)
{
    gr_ctx_t ZZ, Pol, Dop, Pol1, Euler, CC;
    gr_ore_poly_t dop, ldop;
    acb_ode_exponents_t expos;
    acb_t pt, expo;
    acb_poly_t x, f;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;
    ulong flags = 0;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    acb_ode_local_diffop_ctx_init(Pol1, Euler, Dop, ZZ);
    gr_ctx_init_complex_acb(CC, prec);
    gr_ore_poly_init(dop, Dop);
    gr_ore_poly_init(ldop, Euler);
    acb_ode_exponents_init(expos);
    acb_init(pt);
    acb_init(expo);
    acb_poly_init(x);
    acb_poly_init(f);
    acb_mat_init(mat, 2, 2);

    acb_set_ui(pt, 1);
    acb_div_ui(pt, pt, 3, prec);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop,
            "(49*x^2 - 21*x)*D^2 - (49*x^2 + 12)*D + 21*x + 12", Dop));
    status |= acb_ode_local_diffop(ldop, &flags, dop, Dop, Euler);
    status |= acb_ode_exponents(expos, ldop, Euler, CC);
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, expos, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "lclm_xpow");

    acb_poly_set_coeff_si(x, 1, 1);
    acb_poly_set_coeff_acb(x, 0, pt);
    for (slong j = 0; j < 2; j++)
    {
        if (acb_contains_zero(expos->groups[j].leader))
        {
            acb_poly_exp_series(f, x, 2, prec);
        }
        else
        {
            acb_set_ui(expo, 3);
            acb_div_ui(expo, expo, 7, prec);
            acb_poly_pow_acb_series(f, x, expo, 2, prec);
        }
        for (slong i = 0; i < 2; i++)
            check_acb(acb_mat_entry(mat, i, j), f->coeffs + i, prec - 32,
                      "lclm_xpow");
    }

    acb_mat_clear(mat);
    acb_poly_clear(f);
    acb_poly_clear(x);
    acb_clear(expo);
    acb_clear(pt);
    acb_ode_exponents_clear(expos);
    gr_ore_poly_clear(ldop, Euler);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(CC);
    gr_ctx_clear(Euler);
    gr_ctx_clear(Pol1);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* sage: ((x*Dx)^2-2-x).numerical_transition_matrix([0,i], 1e-5)
   [ [-1.4237...] + [0.0706...]*I [-0.7959...] + [0.6169...]*I]
   [ [1.4514...] + [-0.1681...]*I  [0.6742...] + [1.2971...]*I]
   Exponents -sqrt(2), sqrt(2); the basis consists of the
   x^a*0F1(; 1 + 2*a; x) for a = -sqrt(2), sqrt(2).
   [ore_algebra] */
static void
test_hyp0f1_sqrt2(void)
{
    gr_ctx_t ZZ, Pol, Dop, Pol1, Euler, CC;
    gr_ore_poly_t dop, ldop;
    acb_ode_exponents_t expos;
    acb_t pt, a, b;
    acb_ptr f;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;
    ulong flags = 0;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    acb_ode_local_diffop_ctx_init(Pol1, Euler, Dop, ZZ);
    gr_ctx_init_complex_acb(CC, prec);
    gr_ore_poly_init(dop, Dop);
    gr_ore_poly_init(ldop, Euler);
    acb_ode_exponents_init(expos);
    acb_init(pt);
    acb_init(a);
    acb_init(b);
    f = _acb_vec_init(2);
    acb_mat_init(mat, 2, 2);

    acb_onei(pt);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "x^2*D^2 + x*D - x - 2", Dop));
    status |= acb_ode_local_diffop(ldop, &flags, dop, Dop, Euler);
    status |= acb_ode_exponents(expos, ldop, Euler, CC);
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, expos, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "hyp0f1_sqrt2");

    for (slong j = 0; j < 2; j++)
    {
        acb_set_si(a, 2);
        acb_sqrt(a, a, prec);
        if (arb_is_negative(acb_realref(expos->groups[j].leader)))
            acb_neg(a, a);
        acb_mul_2exp_si(b, a, 1);
        acb_add_si(b, b, 1, prec);
        xpow_pfq_jet(f, a, NULL, 0, b, 1, pt, 2, prec);
        for (slong i = 0; i < 2; i++)
            check_acb(acb_mat_entry(mat, i, j), f + i, prec - 32,
                      "hyp0f1_sqrt2");
    }

    acb_mat_clear(mat);
    _acb_vec_clear(f, 2);
    acb_clear(b);
    acb_clear(a);
    acb_clear(pt);
    acb_ode_exponents_clear(expos);
    gr_ore_poly_clear(ldop, Euler);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(CC);
    gr_ctx_clear(Euler);
    gr_ctx_clear(Pol1);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* sage: ((x*Dx)^3-2-x).numerical_transition_matrix([0,1], 1e-5)
   [[0.75335...] + [-0.09777...]*I  [0.75335...] + [0.09777...]*I  [1.1080...] + [+/- ...]*I]
   [[-0.78644...] + [-0.8794...]*I [-0.78644...] + [0.8794...]*I  [1.5074...] + [+/- ...]*I]
   [ [0.10146...] + [1.2165...]*I  [0.10146...] + [-1.2165...]*I  [0.32506...] + [+/- ...]*I]
   The exponents are the cube roots a_k of 2, and the basis elements are
   the x^a_k*0F2(; 1 + a_k - a_l, 1 + a_k - a_m; x), {k, l, m} = {0, 1, 2}.
   [ore_algebra] */
static void
test_hyp0f2_cbrt2(void)
{
    gr_ctx_t ZZ, Pol, Dop, Pol1, Euler, CC;
    gr_ore_poly_t dop, ldop;
    acb_ode_exponents_t expos;
    acb_t pt, c;
    acb_ptr expo, b;
    acb_ptr f;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;
    ulong flags = 0;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    acb_ode_local_diffop_ctx_init(Pol1, Euler, Dop, ZZ);
    gr_ctx_init_complex_acb(CC, prec);
    gr_ore_poly_init(dop, Dop);
    gr_ore_poly_init(ldop, Euler);
    acb_ode_exponents_init(expos);
    acb_init(pt);
    acb_init(c);
    expo = _acb_vec_init(3);
    b = _acb_vec_init(2);
    f = _acb_vec_init(3);
    acb_mat_init(mat, 3, 3);

    acb_one(pt);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop,
            "x^3*D^3 + 3*x^2*D^2 + x*D - x - 2", Dop));
    status |= acb_ode_local_diffop(ldop, &flags, dop, Dop, Euler);
    status |= acb_ode_exponents(expos, ldop, Euler, CC);
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, expos, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "hyp0f2_cbrt2");

    /* Im(a_0) < 0 < Im(a_1), a_2 real */
    set_acb_str(c, "2^(1/3)", prec);
    set_acb_str(expo + 0, "exp(-2*pi*I/3)", prec);
    acb_mul(expo + 0, expo + 0, c, prec);
    acb_conj(expo + 1, expo + 0);
    acb_set(expo + 2, c);

    for (slong j = 0; j < 3; j++)
    {
        const arb_struct * im = acb_imagref(expos->groups[j].leader);
        slong k = arb_is_negative(im) ? 0 : arb_is_positive(im) ? 1 : 2;

        for (slong l = 0, i = 0; l < 3; l++)
        {
            if (l == k)
                continue;
            acb_sub(b + i, expo + k, expo + l, prec);
            acb_add_si(b + i, b + i, 1, prec);
            i++;
        }
        xpow_pfq_jet(f, expo + k, NULL, 0, b, 2, pt, 3, prec);
        for (slong i = 0; i < 3; i++)
            check_acb(acb_mat_entry(mat, i, j), f + i, prec - 32,
                      "hyp0f2_cbrt2");
    }

    acb_mat_clear(mat);
    _acb_vec_clear(f, 3);
    _acb_vec_clear(b, 2);
    _acb_vec_clear(expo, 3);
    acb_clear(c);
    acb_clear(pt);
    acb_ode_exponents_clear(expos);
    gr_ore_poly_clear(ldop, Euler);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(CC);
    gr_ctx_clear(Euler);
    gr_ctx_clear(Pol1);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* E_2(z) = cosh(sqrt(z)), as the element of exponent 0 of the local basis at
   0 of 4*x^2*Dx^2 + 2*x*Dx - x (exponents 0, 1/2) */
static int
my_mittag_leffler_2(acb_t res, const acb_t z, slong prec)
{
    gr_ctx_t ZZ, Pol, Dop, Pol1, Euler, CC;
    gr_ore_poly_t dop, ldop;
    acb_ode_exponents_t expos;
    acb_mat_t mat;
    int status = GR_SUCCESS;
    ulong flags = 0;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    acb_ode_local_diffop_ctx_init(Pol1, Euler, Dop, ZZ);
    gr_ctx_init_complex_acb(CC, prec);
    gr_ore_poly_init(dop, Dop);
    gr_ore_poly_init(ldop, Euler);
    acb_ode_exponents_init(expos);
    acb_mat_init(mat, 1, 2);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "4*x^2*D^2 + 2*x*D - x", Dop));
    status |= acb_ode_local_diffop(ldop, &flags, dop, Dop, Euler);
    status |= acb_ode_exponents(expos, ldop, Euler, CC);
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, expos, NULL, z,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    if (status == GR_SUCCESS)
    {
        slong j = acb_contains_zero(expos->groups[0].leader) ? 0 : 1;
        acb_swap(res, acb_mat_entry(mat, 0, j));
    }

    acb_mat_clear(mat);
    acb_ode_exponents_clear(expos);
    gr_ore_poly_clear(ldop, Euler);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(CC);
    gr_ctx_clear(Euler);
    gr_ctx_clear(Pol1);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);

    return status;
}

static void
ref_mittag_leffler_2(acb_t res, const acb_t z, slong prec)
{
    acb_sqrt(res, z, prec);
    acb_cosh(res, res, prec);
}

/* sage: e2 = stdfun.mittag_leffler_e(2)
   sage: e2.dop(cosh(sqrt(x))).simplify_full()
   0
   sage: e2.dop.numerical_solution(e2.ini, [e2.pt, 2])
   [2.17818355660857...]
   with e2.dop = 4*(x*Dx)^2 - 2*x*Dx - x (up to a constant factor).
   [ore_algebra] */
static void
test_mittag_leffler_2(void)
{
    acb_t z, res, ref;
    slong prec = 128;

    acb_init(z);
    acb_init(res);
    acb_init(ref);

    acb_set_si(z, 2);
    check_status(my_mittag_leffler_2(res, z, prec), "mittag_leffler_2");
    ref_mittag_leffler_2(ref, z, prec);
    check_acb(res, ref, prec - 32, "mittag_leffler_2");

    acb_clear(ref);
    acb_clear(res);
    acb_clear(z);
}

/* E_{1/2}(z) as a combination of the local basis at 0 of Dx^2 - 2*x*Dx - 2 */
static int
my_mittag_leffler_half(acb_t res, const acb_t z, slong prec)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_mat_t mat;
    acb_t c;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_mat_init(mat, 1, 2);
    acb_init(c);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "D^2 - 2*x*D - 2", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, z,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    set_acb_str(c, "2/sqrt(pi)", prec);
    acb_mul(res, acb_mat_entry(mat, 0, 1), c, prec);
    acb_add(res, res, acb_mat_entry(mat, 0, 0), prec);

    acb_clear(c);
    acb_mat_clear(mat);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);

    return status;
}

/* E_{1/2}(z) = exp(z^2)*(1 + erf(z)) */
static void
ref_mittag_leffler_half(acb_t res, const acb_t z, slong prec)
{
    acb_t c;

    acb_init(c);

    acb_hypgeom_erf(res, z, prec);
    acb_add_si(res, res, 1, prec);
    acb_sqr(c, z, prec);
    acb_exp(c, c, prec);
    acb_mul(res, res, c, prec);

    acb_clear(c);
}

/* sage: eh = stdfun.mittag_leffler_e(1/2)
   sage: eh.dop(exp(x^2)*(1+erf(x))).simplify_full()
   0
   sage: eh.dop.numerical_solution(eh.ini, [eh.pt, 2])
   [108.94090438997...]
   with eh.dop = (x*Dx - 1)*(x*Dx - 2*x^2) = x^2*(Dx^2 - 2*x*Dx - 2) and
   eh.ini = [1, 2/sqrt(pi)].
   [ore_algebra] */
static void
test_mittag_leffler_half(void)
{
    acb_t z, res, ref;
    slong prec = 128;

    acb_init(z);
    acb_init(res);
    acb_init(ref);

    acb_set_si(z, 2);
    check_status(my_mittag_leffler_half(res, z, prec),
                 "mittag_leffler_half");
    ref_mittag_leffler_half(ref, z, prec);
    check_acb(res, ref, prec - 32, "mittag_leffler_half");

    acb_clear(ref);
    acb_clear(res);
    acb_clear(z);
}

/* log(z), as the coefficient of log(x) in the local basis log(x), 1 at 0 of
   Dx*x*Dx = x*Dx^2 + Dx */
static int
my_log(acb_t res, const acb_t z, slong prec)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_mat_t mat;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_mat_init(mat, 1, 2);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "x*D^2 + D", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, z,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    acb_swap(res, acb_mat_entry(mat, 0, 0));

    acb_mat_clear(mat);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);

    return status;
}

/* sage: stdfun.log.dop.numerical_solution(stdfun.log.ini,
   ....:                                   [stdfun.log.pt, 2])
   [0.69314718055994...]
   with stdfun.log.dop = Dx*x*Dx. Here the solution is expanded at 0 rather
   than at stdfun.log.pt = 1.
   [ore_algebra] */
static void
test_log(void)
{
    acb_t z, res, ref;
    slong prec = 128;

    acb_init(z);
    acb_init(res);
    acb_init(ref);

    acb_set_si(z, 2);
    check_status(my_log(res, z, prec), "log");
    acb_log(ref, z, prec);
    check_acb(res, ref, prec - 16, "log");

    acb_clear(ref);
    acb_clear(res);
    acb_clear(z);
}

/* sage: a = QQbar(2*exp(I*2*pi/3))
   sage: b = a + a.conjugate()
   sage: (Dx*x*Dx).numerical_transition_matrix([0, b])[0,0]
   [0.6931471805599453...] + [3.141592653589793...]*I
   Evaluation on the branch cut (b = -2).
   [ore_algebra] */
static void
test_log_on_cut(void)
{
    acb_t z, res, ref;
    slong prec = 128;

    acb_init(z);
    acb_init(res);
    acb_init(ref);

    acb_set_si(z, -2);
    check_status(my_log(res, z, prec), "log_on_cut");
    set_acb_str(ref, "log(2) + pi*I", prec);
    check_acb(res, ref, prec - 16, "log_on_cut");

    acb_clear(ref);
    acb_clear(res);
    acb_clear(z);
}

/* sage: (Dx*x*Dx).numerical_transition_matrix([0, b - I/10^1000])[0,0]
   [0.6931471805599453...] + [-3.141592653589793...]*I
   (continued from test_log_on_cut, with 2^-1000 in place of 10^-1000)
   [ore_algebra] */
static void
test_log_below_cut(void)
{
    acb_t z, res, ref;
    slong prec = 128;

    acb_init(z);
    acb_init(res);
    acb_init(ref);

    acb_set_si(z, -2);
    arb_set_si(acb_imagref(z), -1);
    arb_mul_2exp_si(acb_imagref(z), acb_imagref(z), -1000);
    check_status(my_log(res, z, prec), "log_below_cut");
    set_acb_str(ref, "log(2) - pi*I", prec);
    check_acb(res, ref, prec - 16, "log_below_cut");

    acb_clear(ref);
    acb_clear(res);
    acb_clear(z);
}

/****************************************************************************
 * Other expansion points, analytic continuation
 ****************************************************************************/

/* The operator (x^2 + 1)*Dx^2 + 2*x*Dx of test_arctan shifted by 1,
   evaluated at 1/2. The basis is 1, 2*(arctan(1 + x) - pi/4). */
static void
test_arctan_shifted(void)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    acb_t pt, c;
    acb_poly_t x, f;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    acb_init(pt);
    acb_init(c);
    acb_poly_init(x);
    acb_poly_init(f);
    acb_mat_init(mat, 2, 2);

    acb_set_d(pt, 0.5);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop,
            "(x^2 + 2*x + 2)*D^2 + (2*x + 2)*D", Dop));
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "arctan_shifted");

    acb_poly_set_coeff_si(x, 1, 1);
    acb_add_si(c, pt, 1, prec);
    acb_poly_set_coeff_acb(x, 0, c);
    acb_poly_atan_series(f, x, 2, prec);
    acb_const_pi(c, prec);
    acb_mul_2exp_si(c, c, -2);
    acb_sub(f->coeffs, f->coeffs, c, prec);
    acb_poly_scalar_mul_2exp_si(f, f, 1);
    for (slong i = 0; i < 2; i++)
    {
        acb_set_si(c, i == 0);
        check_acb(acb_mat_entry(mat, i, 0), c, prec - 16, "arctan_shifted");
        check_acb(acb_mat_entry(mat, i, 1), f->coeffs + i, prec - 16,
                  "arctan_shifted");
    }

    acb_mat_clear(mat);
    acb_poly_clear(f);
    acb_poly_clear(x);
    acb_clear(c);
    acb_clear(pt);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* res(x) = dop(x + a), i.e., dop rewritten in terms of the local variable
   x - a */
static void
shift_dop(gr_ore_poly_t res, const gr_ore_poly_t dop, gr_srcptr a,
          gr_ore_poly_ctx_t Dop)
{
    gr_ctx_struct * Cst;

    GR_MUST_SUCCEED(gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Cst, NULL, Dop));

    GR_MUST_SUCCEED(gr_ore_poly_set(res, dop, Dop));
    for (slong i = 0; i < res->length; i++)
    {
        gr_poly_struct * c = gr_ore_poly_coeff_ptr(res, i, Dop);
        GR_MUST_SUCCEED(gr_poly_taylor_shift(c, c, a, Cst));
    }
}

/* sage: r = QQbar(I*sqrt(2))
   sage: b = a + a.conjugate() + r
   sage: (Dx*(x^2+2)*Dx).numerical_transition_matrix([r, b])[0,0]
   [0.4904146265058631...] + [2.526112944919405...]*I
   (with a as in test_log_on_cut, so that b - r = -2.) Expansion at the
   singular point i*sqrt(2). In terms of t = x - i*sqrt(2), the basis is
   log(t) - log(1 + t/(2*i*sqrt(2))), 1.
   [ore_algebra] */
static void
test_log_algebraic_singularity(void)
{
    gr_ctx_t QQbar, Pol, Dop;
    gr_ore_poly_t dop;
    gr_ptr a;
    acb_t t, c;
    acb_poly_t x, f, g;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_complex_qqbar(QQbar);
    gr_ctx_init_gr_poly(Pol, QQbar);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    a = gr_heap_init(QQbar);
    acb_init(t);
    acb_init(c);
    acb_poly_init(x);
    acb_poly_init(f);
    acb_poly_init(g);
    acb_mat_init(mat, 2, 2);

    acb_set_si(t, -2);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "(x^2 + 2)*D^2 + 2*x*D", Dop));
    GR_MUST_SUCCEED(gr_set_str(a, "I*sqrt(2)", QQbar));
    shift_dop(dop, dop, a, Dop);
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, t,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);
    check_status(status, "log_algebraic_singularity");

    acb_poly_set_coeff_si(x, 1, 1);
    acb_poly_set_coeff_acb(x, 0, t);
    acb_poly_log_series(f, x, 2, prec);
    set_acb_str(c, "1/(2*I*sqrt(2))", prec);
    acb_poly_scalar_mul(g, x, c, prec);
    acb_poly_log1p_series(g, g, 2, prec);
    acb_poly_sub(f, f, g, prec);
    for (slong i = 0; i < 2; i++)
    {
        check_acb(acb_mat_entry(mat, i, 0), f->coeffs + i, prec - 32,
                  "log_algebraic_singularity");
        acb_set_si(c, i == 0);
        check_acb(acb_mat_entry(mat, i, 1), c, prec - 32,
                  "log_algebraic_singularity");
    }

    acb_mat_clear(mat);
    acb_poly_clear(g);
    acb_poly_clear(f);
    acb_poly_clear(x);
    acb_clear(c);
    acb_clear(t);
    gr_heap_clear(a, QQbar);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(QQbar);
}

/* Transition matrix of dop_str, an operator in QQbar[x]<D>, along the
   polygonal path with vertices path[0], ..., path[len-1] (elements of QQbar
   given as strings). The vertices must be ordinary points, and each step
   must stay within the disk of convergence of the local solutions at its
   origin. */
static int
transition_matrix_path(acb_mat_t mat, const char * dop_str,
                       const char ** path, slong len, slong prec)
{
    gr_ctx_t QQbar, CC, Pol, Dop;
    gr_ore_poly_t dop, sdop;
    gr_ptr a, b;
    acb_t pt;
    acb_mat_t step;
    int status = GR_SUCCESS;

    gr_ctx_init_complex_qqbar(QQbar);
    gr_ctx_init_complex_acb(CC, prec);
    gr_ctx_init_gr_poly(Pol, QQbar);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    gr_ore_poly_init(sdop, Dop);
    a = gr_heap_init(QQbar);
    b = gr_heap_init(QQbar);
    acb_init(pt);
    acb_mat_init(step, acb_mat_nrows(mat), acb_mat_ncols(mat));

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, dop_str, Dop));
    GR_MUST_SUCCEED(gr_set_str(a, path[0], QQbar));
    acb_mat_one(mat);

    for (slong k = 1; k < len && status == GR_SUCCESS; k++)
    {
        GR_MUST_SUCCEED(gr_set_str(b, path[k], QQbar));
        shift_dop(sdop, dop, a, Dop);
        GR_MUST_SUCCEED(gr_sub(a, b, a, QQbar));
        GR_MUST_SUCCEED(gr_set_other(pt, a, QQbar, CC));
        status |= acb_ode_fundamental_matrix(step, sdop, Dop, NULL, NULL, pt,
                                             ACB_ODE_BASIS_ECHELON, 0, prec);
        acb_mat_mul(mat, step, mat, prec);
        gr_swap(a, b, QQbar);
    }

    acb_mat_clear(step);
    acb_clear(pt);
    gr_heap_clear(b, QQbar);
    gr_heap_clear(a, QQbar);
    gr_ore_poly_clear(sdop, Dop);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(CC);
    gr_ctx_clear(QQbar);

    return status;
}

/* sage: dop = (x^2 + 1)*Dx^2 + 2*x*Dx
   sage: dop.numerical_transition_matrix([0, 1], 1e-20)
   [  [1.00...] [0.7853981633974483096...]]
   [  [+/- ...] [0.5000000000000000000...]]
   The basis is 1, arctan(x). The evaluation point lies on the boundary of
   the disk of convergence at 0, hence the intermediate point 1/2.
   [ore_algebra] */
static void
test_arctan(void)
{
    const char * path[] = { "0", "1/2", "1" };
    acb_t pt, c;
    acb_poly_t x, f;
    acb_mat_t mat;
    slong prec = 128;

    acb_init(pt);
    acb_init(c);
    acb_poly_init(x);
    acb_poly_init(f);
    acb_mat_init(mat, 2, 2);

    check_status(transition_matrix_path(mat, "(x^2 + 1)*D^2 + 2*x*D",
                                        path, 3, prec), "arctan");

    acb_one(pt);
    acb_poly_set_coeff_si(x, 1, 1);
    acb_poly_set_coeff_acb(x, 0, pt);
    acb_poly_atan_series(f, x, 2, prec);
    for (slong i = 0; i < 2; i++)
    {
        acb_set_si(c, i == 0);
        check_acb(acb_mat_entry(mat, i, 0), c, prec - 16, "arctan");
        check_acb(acb_mat_entry(mat, i, 1), f->coeffs + i, prec - 16,
                  "arctan");
    }

    acb_mat_clear(mat);
    acb_poly_clear(f);
    acb_poly_clear(x);
    acb_clear(c);
    acb_clear(pt);
}

/* sage: dop.numerical_transition_matrix([0, 1+i], 1e-10)
   [ [1.00...] + [+/- ...]*I  [1.017221967...] + [0.4023594781...]*I]
   [ [+/- ...] + [+/- ...]*I  [0.200000000...] + [-0.400000000...]*I]
   (continued from test_arctan)
   [ore_algebra] */
static void
test_arctan_complex(void)
{
    const char * path[] = { "0", "1/2", "1 + I/2", "1 + I" };
    acb_t pt, c;
    acb_poly_t x, f;
    acb_mat_t mat;
    slong prec = 128;

    acb_init(pt);
    acb_init(c);
    acb_poly_init(x);
    acb_poly_init(f);
    acb_mat_init(mat, 2, 2);

    check_status(transition_matrix_path(mat, "(x^2 + 1)*D^2 + 2*x*D",
                                        path, 4, prec), "arctan_complex");

    set_acb_str(pt, "1 + I", prec);
    acb_poly_set_coeff_si(x, 1, 1);
    acb_poly_set_coeff_acb(x, 0, pt);
    acb_poly_atan_series(f, x, 2, prec);
    for (slong i = 0; i < 2; i++)
    {
        acb_set_si(c, i == 0);
        check_acb(acb_mat_entry(mat, i, 0), c, prec - 16, "arctan_complex");
        check_acb(acb_mat_entry(mat, i, 1), f->coeffs + i, prec - 16,
                  "arctan_complex");
    }

    acb_mat_clear(mat);
    acb_poly_clear(f);
    acb_poly_clear(x);
    acb_clear(c);
    acb_clear(pt);
}

/* sage: dop = (x+1)*(x^2+1)*Dx^3-(x-1)*(x^2-3)*Dx^2-2*(x^2+2*x-1)*Dx
   sage: dop.numerical_transition_matrix( [0, 1+i], 1e-12)
   [ [1.00...] + [+/- ...]*I  [1.017221967...] + [0.402359478...]*I  [-1.097056...] + [3.76999161...]*I]
   [ [+/- ...] + [+/- ...]*I  [0.200000000...] + [-0.40000000...]*I  [2.5373878...] + [5.37471057...]*I]
   [ [+/- ...] + [+/- ...]*I  [-0.04000000...] + [0.280000000...]*I  [1.5486939...] + [1.72735528...]*I]
   The basis is 1, arctan(x), 2*(exp(x) - 1 - arctan(x)).
   [ore_algebra] */
static void
test_exp_arctan(void)
{
    const char * path[] = { "0", "1/2", "1", "1 + I/2", "1 + I" };
    acb_t pt, c;
    acb_poly_t x, f, g;
    acb_mat_t mat;
    slong prec = 128;

    acb_init(pt);
    acb_init(c);
    acb_poly_init(x);
    acb_poly_init(f);
    acb_poly_init(g);
    acb_mat_init(mat, 3, 3);

    check_status(transition_matrix_path(mat,
            "(x + 1)*(x^2 + 1)*D^3 - (x - 1)*(x^2 - 3)*D^2"
            " - 2*(x^2 + 2*x - 1)*D", path, 5, prec), "exp_arctan");

    set_acb_str(pt, "1 + I", prec);
    acb_poly_set_coeff_si(x, 1, 1);
    acb_poly_set_coeff_acb(x, 0, pt);
    acb_poly_atan_series(f, x, 3, prec);
    acb_poly_exp_series(g, x, 3, prec);
    acb_poly_sub(g, g, f, prec);
    acb_poly_add_si(g, g, -1, prec);
    acb_poly_scalar_mul_2exp_si(g, g, 1);
    for (slong i = 0; i < 3; i++)
    {
        acb_set_si(c, i == 0);
        check_acb(acb_mat_entry(mat, i, 0), c, prec - 16, "exp_arctan");
        check_acb(acb_mat_entry(mat, i, 1), f->coeffs + i, prec - 16,
                  "exp_arctan");
        check_acb(acb_mat_entry(mat, i, 2), g->coeffs + i, prec - 16,
                  "exp_arctan");
    }

    acb_mat_clear(mat);
    acb_poly_clear(g);
    acb_poly_clear(f);
    acb_poly_clear(x);
    acb_clear(c);
    acb_clear(pt);
}

/* sage: (Dx*x*Dx).numerical_transition_matrix([1, i, -1, -i, 1])
   [ [1.0000...] + [+/- ...]*I [+/- ...] + [6.2831853071795...]*I]
   [   [+/- ...] + [+/- ...]*I          [1.0000...] + [+/- ...]*I]
   Monodromy of log(x) around the origin, with additional vertices so that
   each step stays inside the disk of convergence.
   [ore_algebra] */
static void
test_log_monodromy(void)
{
    const char * path[] = {
        "1", "3/4 + 3/4*I", "I", "-3/4 + 3/4*I", "-1",
        "-3/4 - 3/4*I", "-I", "3/4 - 3/4*I", "1",
    };
    acb_t ref;
    acb_mat_t mat;
    slong prec = 128;

    acb_init(ref);
    acb_mat_init(mat, 2, 2);

    check_status(transition_matrix_path(mat, "x*D^2 + D", path, 9, prec),
                 "log_monodromy");

    for (slong i = 0; i < 2; i++)
    {
        for (slong j = 0; j < 2; j++)
        {
            if (i == 0 && j == 1)
                set_acb_str(ref, "2*pi*I", prec);
            else
                acb_set_si(ref, i == j);
            check_acb(acb_mat_entry(mat, i, j), ref, prec - 32,
                      "log_monodromy");
        }
    }

    acb_mat_clear(mat);
    acb_clear(ref);
}

/* sage: ode = D + 1/4/(t - 1/2)
   sage: ode.numerical_transition_matrix([0,1+I,1], 1e-110,
   ....:                                 algorithm='binsplit')
   [[0.707...2078...] + [0.707...]*I]
   The solution is (1 - 2*t)^(-1/4); continuing it from 0 to 1 through the
   upper half plane yields exp(i*pi/4). The operator is scaled to
   (4*x - 2)*Dx + 1 and the path subdivided around the singular point 1/2.
   [ore_algebra] */
static void
test_quartic_branch(void)
{
    const char * path[] = {
        "0", "I/4", "I/2", "1/2 + 3/4*I", "1 + I/2", "1",
    };
    acb_t ref;
    acb_mat_t mat;
    slong prec = 384;

    acb_init(ref);
    acb_mat_init(mat, 1, 1);

    check_status(transition_matrix_path(mat, "(4*x - 2)*D + 1", path, 6,
                                        prec), "quartic_branch");

    set_acb_str(ref, "exp(pi*I/4)", prec);
    check_acb(acb_mat_entry(mat, 0, 0), ref, prec - 96, "quartic_branch");

    acb_mat_clear(mat);
    acb_clear(ref);
}

/* sage: dop = (z+1)*(3*z^2-z+2)*Dz^3 + (5*z^3+4*z^2+2*z+4)*Dz^2 \
   ....:       + (z+1)*Dz + (4*z^3+2*z^2+5)
   sage: path = [0,-2/5+3/5*i,-2/5+i,-1/5+7/5*i]
   sage: dop.numerical_solution([0,i,0], path, 1e-150)
   [-1.5598481440603221187326507993405933893413346644879595004537063375459901302359572361012065551669069...] +
   [-0.7107764943512671843673286878693314397759047479618104045777076954591551406949345143368742955333566...]*I
   The path is subdivided so that each step has a ratio below 1/2 and, from
   -2/5 + 4/5*i on, below 1/5. (At ratios of about 0.35, the steps starting
   from -2/5 + 4/5*i, -2/5 + i, and -3/10 + 6/5*i lose 25-35% of the working
   precision.) The reference value only has about 320 bits, but we check
   that the result is accurate to the requested 1e-150.
   [ore_algebra] */
static void
test_order3_high_prec(void)
{
    const char * path[] = {
        "0", "-2/15 + 1/5*I", "-4/15 + 2/5*I", "-2/5 + 3/5*I",
        "-2/5 + 4/5*I", "-2/5 + 9/10*I", "-2/5 + I",
        "-7/20 + 11/10*I", "-3/10 + 6/5*I", "-1/4 + 13/10*I", "-1/5 + 7/5*I",
    };
    acb_t ref, res;
    acb_mat_t mat;
    slong prec = 544;

    acb_init(ref);
    acb_init(res);
    acb_mat_init(mat, 3, 3);

    check_status(transition_matrix_path(mat,
            "(x + 1)*(3*x^2 - x + 2)*D^3 + (5*x^3 + 4*x^2 + 2*x + 4)*D^2"
            " + (x + 1)*D + 4*x^3 + 2*x^2 + 5", path, 11, prec),
            "order3_high_prec");

    acb_mul_onei(res, acb_mat_entry(mat, 0, 1));
    set_acb_str(ref,
        "[-1.5598481440603221187326507993405933893413346644879595004537063375"
        "459901302359572361012065551669069 +/- 1e-97]"
        " + [-0.710776494351267184367328687869331439775904747961810404577707"
        "6954591551406949345143368742955333566 +/- 1e-97]*I", prec);
    check_acb(res, ref, 500, "order3_high_prec");

    acb_mat_clear(mat);
    acb_clear(res);
    acb_clear(ref);
}

/* The operator dop4 annihilating the generating series of face-centered
   cubic lattice walks (4d), from the examples of ore_algebra. */
#define FCC_DOP4 \
    "(-1 + x)*x^3*(2 + x)*(3 + x)*(6 + x)*(8 + x)*(4 + 3*x)^2*D^4" \
    " + 2*x^2*(4 + 3*x)*(-3456 - 2304*x + 3676*x^2 + 4920*x^3 + 2079*x^4" \
    " + 356*x^5 + 21*x^6)*D^3" \
    " + 6*x*(-5376 - 5248*x + 11080*x^2 + 25286*x^3 + 19898*x^4" \
    " + 7432*x^5 + 1286*x^6 + 81*x^7)*D^2" \
    " + 12*(-384 + 224*x + 3716*x^2 + 7633*x^3 + 6734*x^4 + 2939*x^5" \
    " + 604*x^6 + 45*x^7)*D" \
    " + 12*x*(256 + 632*x + 702*x^2 + 382*x^3 + 98*x^4 + 9*x^5)"

/* Fundamental matrix at 1/2 of FCC_DOP4 in its local basis at a */
static int
fcc4(acb_mat_t mat, slong a, slong prec)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ore_poly_t dop;
    gr_ptr c;
    acb_t pt;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    c = gr_heap_init(ZZ);
    acb_init(pt);

    acb_set_d(pt, 0.5);
    acb_sub_si(pt, pt, a, prec);
    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, FCC_DOP4, Dop));
    GR_MUST_SUCCEED(gr_set_si(c, a, ZZ));
    shift_dop(dop, dop, c, Dop);
    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);

    acb_clear(pt);
    gr_heap_clear(c, ZZ);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);

    return status;
}

/* sage: fcc.dop4.numerical_solution([0, 1, 0, 0], [1, 1/2],
   ....:                             algorithm="binsplit")
   [0.6757976924611368...] + [-2.987541622123742...]*I
   Expansion at the regular singular point 1 (exponents 0, 1, 1, 2).
   [ore_algebra] */
static void
test_fcc4_singular_base(void)
{
    acb_t ref;
    acb_mat_t mat;
    slong prec = 128;

    acb_init(ref);
    acb_mat_init(mat, 4, 4);

    check_status(fcc4(mat, 1, prec), "fcc4_singular_base");

    set_acb_str(ref, "[0.6757976924611368 +/- 1e-16]"
                     " + [-2.987541622123742 +/- 1e-15]*I", prec);
    check_acb(acb_mat_entry(mat, 0, 1), ref, prec - 48, "fcc4_singular_base");

    acb_mat_clear(mat);
    acb_clear(ref);
}

/* sage: fcc.dop4.numerical_solution([0, 0, 0, 1], [0, 1], 1e-60)
   [1.10584379792120476018299547088585107443954623663875285836...]
   + [+/- ...]*I
   Connection between the regular singular points 0 and 1: solve
   M1*y = M0*[0, 0, 0, 1], where M0, M1 are the fundamental matrices of the
   local bases at 0 and 1, both evaluated at 1/2. The value at 1 is the
   coordinate y[0] on the basis element of leading monomial 1.
   [ore_algebra] */
static void
test_fcc4_connection(void)
{
    acb_t ref;
    acb_mat_t mat0, mat1, rhs, y;
    slong prec = 192;

    acb_init(ref);
    acb_mat_init(mat0, 4, 4);
    acb_mat_init(mat1, 4, 4);
    acb_mat_init(rhs, 4, 1);
    acb_mat_init(y, 4, 1);

    check_status(fcc4(mat0, 0, prec), "fcc4_connection");
    check_status(fcc4(mat1, 1, prec), "fcc4_connection");

    for (slong i = 0; i < 4; i++)
        acb_set(acb_mat_entry(rhs, i, 0), acb_mat_entry(mat0, i, 3));
    if (!acb_mat_solve(y, mat1, rhs, prec))
    {
        flint_printf("FAIL: fcc4_connection (singular matrix)\n");
        flint_abort();
    }

    set_acb_str(ref,
        "[1.10584379792120476018299547088585107443954623663875285836"
        " +/- 1e-56] + [+/- 1e-48]*I", prec);
    check_acb(acb_mat_entry(y, 0, 0), ref, prec - 96, "fcc4_connection");

    acb_mat_clear(y);
    acb_mat_clear(rhs);
    acb_mat_clear(mat1);
    acb_mat_clear(mat0);
    acb_clear(ref);
}

/* A random operator with prescribed singularities and exponents (from
   Chyzak, Goyer, Mezzarobba, 2022) that defines transition matrices with
   large entries, from the doctests of ore_algebra. Its leading coefficient
   is c*(x - 75/17)^3*(x - 4)^3. */
#define CGM_DOP \
    "(21829048560*x^6 - 550862460720*x^5 + 5789380916880*x^4" \
    " - 32434807101840*x^3 + 102165545592000*x^2 - 171548863200000*x" \
    " + 119964240000000)*D^3 + (-506769209364*x^5 + 10798395588024*x^4" \
    " - 91979496877620*x^3 + 391481933812992*x^2 - 832564302926400*x" \
    " + 707776470720000)*D^2 + (251156588660*x^4 - 10005924091768*x^3" \
    " + 99364591554260*x^2 - 379332696633792*x + 503580376843200)*D" \
    " + 38987126730300*x^3 - 113173203800613*x^2 - 1127681839519020*x" \
    " + 3827710032458400"

/* Fundamental matrix at 143/34 of CGM_DOP in its local basis at a, with the
   exponent groups sorted by increasing real part as in ore_algebra */
static int
large_entries(acb_mat_t mat, const char * a_str, slong prec)
{
    gr_ctx_t QQbar, CC, Pol, Dop, Pol1, Euler;
    gr_ore_poly_t dop, ldop;
    acb_ode_exponents_t expos;
    acb_ode_group_struct * grp;
    gr_ptr a, b;
    acb_t pt;
    int status = GR_SUCCESS;
    ulong flags = 0;

    gr_ctx_init_complex_qqbar(QQbar);
    gr_ctx_init_complex_acb(CC, prec);
    gr_ctx_init_gr_poly(Pol, QQbar);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    acb_ode_local_diffop_ctx_init(Pol1, Euler, Dop, QQbar);
    gr_ore_poly_init(dop, Dop);
    gr_ore_poly_init(ldop, Euler);
    acb_ode_exponents_init(expos);
    a = gr_heap_init(QQbar);
    b = gr_heap_init(QQbar);
    acb_init(pt);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, CGM_DOP, Dop));
    GR_MUST_SUCCEED(gr_set_str(a, a_str, QQbar));
    GR_MUST_SUCCEED(gr_set_str(b, "143/34", QQbar));
    GR_MUST_SUCCEED(gr_sub(b, b, a, QQbar));
    GR_MUST_SUCCEED(gr_set_other(pt, b, QQbar, CC));
    shift_dop(dop, dop, a, Dop);
    status |= acb_ode_local_diffop(ldop, &flags, dop, Dop, Euler);
    status |= acb_ode_exponents(expos, ldop, Euler, CC);

    grp = expos->groups;
    for (slong i = 1; i < expos->ngroups; i++)
        for (slong j = i; j > 0 && arb_lt(acb_realref(grp[j].leader),
                                          acb_realref(grp[j - 1].leader)); j--)
            FLINT_SWAP(acb_ode_group_struct, grp[j], grp[j - 1]);

    status |= acb_ode_fundamental_matrix(mat, dop, Dop, expos, NULL, pt,
                                         ACB_ODE_BASIS_ECHELON, 0, prec);

    acb_clear(pt);
    gr_heap_clear(b, QQbar);
    gr_heap_clear(a, QQbar);
    acb_ode_exponents_clear(expos);
    gr_ore_poly_clear(ldop, Euler);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Euler);
    gr_ctx_clear(Pol1);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(CC);
    gr_ctx_clear(QQbar);

    return status;
}

/* sage: dop.numerical_transition_matrix([75/17,4]).trace()
   [-2.048077478517672...e+95...] + [-2.818936813156118...e+95...]*I
   with dop = CGM_DOP. Both 75/17 and 4 are regular singular points. As in
   test_fcc4_connection, we compute the transition matrix as M1^(-1)*M0,
   where M0, M1 are the fundamental matrices of the local bases at 75/17 and
   4 evaluated at the midpoint 143/34. Both evaluations lose about half of
   the working precision (presumably due to the triple roots), and the
   linear solve about 200 more bits, hence the high precision.
   [ore_algebra] */
static void
test_large_entries(void)
{
    acb_t ref, tr;
    acb_mat_t mat0, mat1, mat;
    slong prec = 768;

    acb_init(ref);
    acb_init(tr);
    acb_mat_init(mat0, 3, 3);
    acb_mat_init(mat1, 3, 3);
    acb_mat_init(mat, 3, 3);

    check_status(large_entries(mat0, "75/17", prec), "large_entries");
    check_status(large_entries(mat1, "4", prec), "large_entries");

    if (!acb_mat_solve(mat, mat1, mat0, prec))
    {
        flint_printf("FAIL: large_entries (singular matrix)\n");
        flint_abort();
    }
    acb_mat_trace(tr, mat, prec);

    set_acb_str(ref, "[-2.048077478517672e+95 +/- 1e+80]"
                     " + [-2.818936813156118e+95 +/- 1e+80]*I", prec);
    check_acb(tr, ref, 64, "large_entries");

    acb_mat_clear(mat);
    acb_mat_clear(mat1);
    acb_mat_clear(mat0);
    acb_clear(tr);
    acb_clear(ref);
}

TEST_FUNCTION_START(acb_ode_fundamental_matrix_vs_ref, state)
{
    test_exp();
    test_exp_complex();
    test_cos_sin();
    test_log1p();
    test_log1p_order3();
    test_ordinary_order3();
    test_dawson();
    test_erf();
    test_erf_imag();
    test_erfi();
    test_airy_ai();
    test_airy_ai_jet();
    test_exp_gaussian();
    test_exp_zeta7();
    test_cosh_sqrt_pi();
    test_exp_i_pi();
    test_neher();
    test_neher_3_2();
    test_high_degree();
    test_apparent_singularity();

    test_inverse();
    test_euler_pm_i();
    test_large_exponent();
    test_large_negative_exponent();
    test_large_dispersion_exp();
    test_large_dispersion_inv();
    test_bessel();
    test_bessel_i_k();
    test_singular_hardcoded();
    test_ci();
    test_ci_complex();
    test_log2_lclm();
    test_whittaker();
    test_lclm_xpow();
    test_hyp0f1_sqrt2();
    test_hyp0f2_cbrt2();
    test_mittag_leffler_2();
    test_mittag_leffler_half();
    test_log();
    test_log_on_cut();
    test_log_below_cut();

    test_arctan_shifted();
    test_log_algebraic_singularity();
    test_arctan();
    test_arctan_complex();
    test_exp_arctan();
    test_log_monodromy();
    test_quartic_branch();
    test_order3_high_prec();
    test_fcc4_singular_base();
    test_fcc4_connection();
    test_large_entries();

    TEST_FUNCTION_END(state);
}
