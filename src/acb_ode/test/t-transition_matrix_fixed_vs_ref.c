/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Generated using Claude Code (Fable 5, Opus 5.5). */

/* Tests of acb_ode_transition_matrix_fixed against reference values, mostly
   adapted from doctests of the ore_algebra package for SageMath. When the
   operator and the points of the path are taken over the same ring, we use
   the same context for both, so that the steps are computed exactly. */

#include "test_helpers.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_ode.h"
#include "acb_poly.h"
#include "gr.h"
#include "gr_ore_poly.h"
#include "gr_poly.h"
#include "gr_vec.h"
#include "acb_ode_test_helpers.h"

/* sage: r = QQbar(I*sqrt(2))
   sage: b = a + a.conjugate() + r
   sage: (Dx*(x^2+2)*Dx).numerical_transition_matrix([r, b])[0,0]
   [0.4904146265058631...] + [2.526112944919405...]*I
   sage: (Dx*(x^2+2)*Dx).numerical_transition_matrix([r, b - I/10^100])[0,0]
   [0.4904146265058631...] + [-3.757072362260180...]*I
   (with a as in test_log_on_cut of t-fundamental_matrix_vs_ref.c, so that
   b - r = -2.) In terms of t = x - i*sqrt(2), the local basis at the
   singular point i*sqrt(2) is log(t) - log(1 + t/(2*i*sqrt(2))), 1.
   [ore_algebra] */
static void
test_log_algebraic_singularity(void)
{
    const char * path_str[] = {
        "[I*sqrt(2), I*sqrt(2) - 2]",
        "[I*sqrt(2), I*sqrt(2) - 2 - I/10^100]",
    };
    const char * dx_str[] = { "-2", "-2 - I/10^100" };
    gr_ctx_t ZZ, Pol, Dop, QQbar;
    gr_ore_poly_t dop;
    gr_vec_t path;
    acb_t t, c;
    acb_poly_t x, f, g;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_complex_qqbar(QQbar);
    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    gr_vec_init(path, 0, QQbar);
    acb_init(t);
    acb_init(c);
    acb_poly_init(x);
    acb_poly_init(f);
    acb_poly_init(g);
    acb_mat_init(mat, 2, 2);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "(x^2 + 2)*D^2 + 2*x*D", Dop));

    for (slong k = 0; k < 2; k++)
    {
        GR_MUST_SUCCEED(gr_vec_set_str(path, path_str[k], 1, QQbar));
        status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path, QQbar,
                                                  NULL, NULL, 0, prec);
        check_status(status, "log_algebraic_singularity");

        acb_set_str(t, dx_str[k], prec);
        acb_poly_set_coeff_si(x, 1, 1);
        acb_poly_set_coeff_acb(x, 0, t);
        acb_poly_log_series(f, x, 2, prec);
        acb_set_str(c, "1/(2*I*sqrt(2))", prec);
        acb_poly_scalar_mul(g, x, c, prec);
        acb_poly_log1p_series(g, g, 2, prec);
        acb_poly_sub(f, f, g, prec);
        for (slong i = 0; i < 2; i++)
        {
            check_acb(acb_mat_entry(mat, i, 0), f->coeffs + i, prec - 8,
                      "log_algebraic_singularity");
            acb_set_si(c, i == 0);
            check_acb(acb_mat_entry(mat, i, 1), c, prec - 8,
                      "log_algebraic_singularity");
        }
    }

    acb_mat_clear(mat);
    acb_poly_clear(g);
    acb_poly_clear(f);
    acb_poly_clear(x);
    acb_clear(c);
    acb_clear(t);
    gr_vec_clear(path, QQbar);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
    gr_ctx_clear(QQbar);
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
    gr_ctx_t ZZ, Pol, Dop, QQ;
    gr_ore_poly_t dop;
    gr_vec_t path;
    acb_t c;
    acb_poly_t x, f;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ctx_init_fmpq(QQ);
    gr_ore_poly_init(dop, Dop);
    gr_vec_init(path, 0, QQ);
    acb_init(c);
    acb_poly_init(x);
    acb_poly_init(f);
    acb_mat_init(mat, 2, 2);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "(x^2 + 1)*D^2 + 2*x*D", Dop));
    GR_MUST_SUCCEED(gr_vec_set_str(path, "[0, 1/2, 1]", 1, QQ));
    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path, QQ,
                                              NULL, NULL, 0, prec);
    check_status(status, "arctan");

    acb_poly_set_coeff_si(x, 1, 1);
    acb_poly_set_coeff_si(x, 0, 1);
    acb_poly_atan_series(f, x, 2, prec);
    for (slong i = 0; i < 2; i++)
    {
        acb_set_si(c, i == 0);
        check_acb(acb_mat_entry(mat, i, 0), c, prec - 8, "arctan");
        check_acb(acb_mat_entry(mat, i, 1), f->coeffs + i, prec - 8,
                  "arctan");
    }

    acb_mat_clear(mat);
    acb_poly_clear(f);
    acb_poly_clear(x);
    acb_clear(c);
    gr_vec_clear(path, QQ);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(QQ);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* sage: dop.numerical_transition_matrix([0, 1+i], 1e-10)
   [ [1.00...] + [+/- ...]*I  [1.017221967...] + [0.4023594781...]*I]
   [ [+/- ...] + [+/- ...]*I  [0.200000000...] + [-0.400000000...]*I]
   (continued from test_arctan)
   [ore_algebra] */
static void
test_arctan_complex(void)
{
    gr_ctx_t ZZ, Pol, Dop, CC;
    gr_ore_poly_t dop;
    gr_vec_t path;
    acb_t pt, c;
    acb_poly_t x, f;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ctx_init_complex_acb(CC, prec);
    gr_ore_poly_init(dop, Dop);
    gr_vec_init(path, 0, CC);
    acb_init(pt);
    acb_init(c);
    acb_poly_init(x);
    acb_poly_init(f);
    acb_mat_init(mat, 2, 2);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "(x^2 + 1)*D^2 + 2*x*D", Dop));
    GR_MUST_SUCCEED(gr_vec_set_str(path, "[0, 1/2, 1 + I/2, 1 + I]", 1, CC));
    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path, CC,
                                              NULL, NULL, 0, prec);
    check_status(status, "arctan_complex");

    acb_set_str(pt, "1 + I", prec);
    acb_poly_set_coeff_si(x, 1, 1);
    acb_poly_set_coeff_acb(x, 0, pt);
    acb_poly_atan_series(f, x, 2, prec);
    for (slong i = 0; i < 2; i++)
    {
        acb_set_si(c, i == 0);
        check_acb(acb_mat_entry(mat, i, 0), c, prec - 8, "arctan_complex");
        check_acb(acb_mat_entry(mat, i, 1), f->coeffs + i, prec - 8,
                  "arctan_complex");
    }

    acb_mat_clear(mat);
    acb_poly_clear(f);
    acb_poly_clear(x);
    acb_clear(c);
    acb_clear(pt);
    gr_vec_clear(path, CC);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(CC);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
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
    gr_ctx_t CC, Pol, Dop;
    gr_ore_poly_t dop;
    gr_vec_t path;
    acb_t pt, c;
    acb_poly_t x, f, g;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_complex_acb(CC, prec);
    gr_ctx_init_gr_poly(Pol, CC);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    gr_vec_init(path, 0, CC);
    acb_init(pt);
    acb_init(c);
    acb_poly_init(x);
    acb_poly_init(f);
    acb_poly_init(g);
    acb_mat_init(mat, 3, 3);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop,
            "(x + 1)*(x^2 + 1)*D^3 - (x - 1)*(x^2 - 3)*D^2"
            " - 2*(x^2 + 2*x - 1)*D", Dop));
    GR_MUST_SUCCEED(gr_vec_set_str(path, "[0, 1/2, 1, 1 + I/2, 1 + I]", 1,
                                   CC));
    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path, CC,
                                              NULL, NULL, 0, prec);
    check_status(status, "exp_arctan");

    acb_set_str(pt, "1 + I", prec);
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
        check_acb(acb_mat_entry(mat, i, 0), c, prec - 8, "exp_arctan");
        check_acb(acb_mat_entry(mat, i, 1), f->coeffs + i, prec - 8,
                  "exp_arctan");
        check_acb(acb_mat_entry(mat, i, 2), g->coeffs + i, prec - 8,
                  "exp_arctan");
    }

    acb_mat_clear(mat);
    acb_poly_clear(g);
    acb_poly_clear(f);
    acb_poly_clear(x);
    acb_clear(c);
    acb_clear(pt);
    gr_vec_clear(path, CC);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(CC);
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
    gr_ctx_t ZZ, Pol, Dop, CC;
    gr_ore_poly_t dop;
    gr_vec_t path;
    acb_t ref;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ctx_init_complex_acb(CC, prec);
    gr_ore_poly_init(dop, Dop);
    gr_vec_init(path, 0, CC);
    acb_init(ref);
    acb_mat_init(mat, 2, 2);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "x*D^2 + D", Dop));
    GR_MUST_SUCCEED(gr_vec_set_str(path,
            "[1, 3/4 + 3/4*I, I, -3/4 + 3/4*I, -1,"
            " -3/4 - 3/4*I, -I, 3/4 - 3/4*I, 1]", 1, CC));
    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path, CC,
                                              NULL, NULL, 0, prec);
    check_status(status, "log_monodromy");

    for (slong i = 0; i < 2; i++)
    {
        for (slong j = 0; j < 2; j++)
        {
            if (i == 0 && j == 1)
                acb_set_str(ref, "2*pi*I", prec);
            else
                acb_set_si(ref, i == j);
            check_acb(acb_mat_entry(mat, i, j), ref, prec - 8,
                      "log_monodromy");
        }
    }

    acb_mat_clear(mat);
    acb_clear(ref);
    gr_vec_clear(path, CC);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(CC);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
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
    gr_ctx_t RR, Pol, Dop, CC;
    gr_ore_poly_t dop;
    gr_vec_t path;
    acb_t ref;
    acb_mat_t mat;
    slong prec = 384;
    int status = GR_SUCCESS;

    gr_ctx_init_real_arb(RR, prec);
    gr_ctx_init_gr_poly(Pol, RR);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ctx_init_complex_acb(CC, prec);
    gr_ore_poly_init(dop, Dop);
    gr_vec_init(path, 0, CC);
    acb_init(ref);
    acb_mat_init(mat, 1, 1);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "(4*x - 2)*D + 1", Dop));
    GR_MUST_SUCCEED(gr_vec_set_str(path,
            "[0, I/4, I/2, 1/2 + 3/4*I, 1 + I/2, 1]", 1, CC));
    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path, CC,
                                              NULL, NULL, 0, prec);
    check_status(status, "quartic_branch");

    acb_set_str(ref, "exp(pi*I/4)", prec);
    check_acb(acb_mat_entry(mat, 0, 0), ref, prec - 8, "quartic_branch");

    acb_mat_clear(mat);
    acb_clear(ref);
    gr_vec_clear(path, CC);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(CC);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(RR);
}

/* sage: dop = (z+1)*(3*z^2-z+2)*Dz^3 + (5*z^3+4*z^2+2*z+4)*Dz^2 \
   ....:       + (z+1)*Dz + (4*z^3+2*z^2+5)
   sage: path = [0,-2/5+3/5*i,-2/5+i,-1/5+7/5*i]
   sage: dop.numerical_solution([0,i,0], path, 1e-150)
   [-1.5598481440603221187326507993405933893413346644879595004537063375459901302359572361012065551669069...] +
   [-0.7107764943512671843673286878693314397759047479618104045777076954591551406949345143368742955333566...]*I
   The path is subdivided so that each step has a ratio below 1/2 and, from
   -2/5 + 4/5*i on, below 1/5. The reference value only has about 320 bits,
   but we check that the result is accurate to the requested 1e-150.
   [ore_algebra] */
static void
test_order3_high_prec(void)
{
    gr_ctx_t ZZ, Pol, Dop, CC;
    gr_ore_poly_t dop;
    gr_vec_t path;
    acb_t ref, res;
    acb_mat_t mat;
    slong prec = 512;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ctx_init_complex_acb(CC, prec);
    gr_ore_poly_init(dop, Dop);
    gr_vec_init(path, 0, CC);
    acb_init(ref);
    acb_init(res);
    acb_mat_init(mat, 3, 3);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop,
            "(x + 1)*(3*x^2 - x + 2)*D^3 + (5*x^3 + 4*x^2 + 2*x + 4)*D^2"
            " + (x + 1)*D + 4*x^3 + 2*x^2 + 5", Dop));
    GR_MUST_SUCCEED(gr_vec_set_str(path,
            "[0, -2/15 + 1/5*I, -4/15 + 2/5*I, -2/5 + 3/5*I,"
            " -2/5 + 4/5*I, -2/5 + 9/10*I, -2/5 + I,"
            " -7/20 + 11/10*I, -3/10 + 6/5*I, -1/4 + 13/10*I, -1/5 + 7/5*I]",
            1, CC));
    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path, CC,
                                              NULL, NULL, 0, prec);
    check_status(status, "order3_high_prec");

    acb_mul_onei(res, acb_mat_entry(mat, 0, 1));
    acb_set_str(ref,
        "[-1.5598481440603221187326507993405933893413346644879595004537063375"
        "459901302359572361012065551669069 +/- 1e-97]"
        " + [-0.710776494351267184367328687869331439775904747961810404577707"
        "6954591551406949345143368742955333566 +/- 1e-97]*I", prec);
    check_acb(res, ref, prec - 12, "order3_high_prec");

    acb_mat_clear(mat);
    acb_clear(res);
    acb_clear(ref);
    gr_vec_clear(path, CC);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(CC);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
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

/* sage: fcc.dop4.numerical_solution([0, 1, 0, 0], [1, 1/2],
   ....:                             algorithm="binsplit")
   [0.6757976924611368...] + [-2.987541622123742...]*I
   Expansion at the regular singular point 1 (exponents 0, 1, 1, 2).
   [ore_algebra] */
static void
test_fcc4_singular_base(void)
{
    gr_ctx_t ZZ, Pol, Dop, QQ;
    gr_ore_poly_t dop;
    gr_vec_t path;
    acb_t ref;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ctx_init_fmpq(QQ);
    gr_ore_poly_init(dop, Dop);
    gr_vec_init(path, 0, QQ);
    acb_init(ref);
    acb_mat_init(mat, 4, 4);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, FCC_DOP4, Dop));
    GR_MUST_SUCCEED(gr_vec_set_str(path, "[1, 1/2]", 1, QQ));
    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path, QQ,
                                              NULL, NULL, 0, prec);
    check_status(status, "fcc4_singular_base");

    acb_set_str(ref, "[0.6757976924611368 +/- 1e-16]"
                     " + [-2.987541622123742 +/- 1e-15]*I", prec);
    check_acb(acb_mat_entry(mat, 0, 1), ref, prec - 4, "fcc4_singular_base");

    acb_mat_clear(mat);
    acb_clear(ref);
    gr_vec_clear(path, QQ);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(QQ);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* sage: fcc.dop4.numerical_solution([0, 0, 0, 1], [0, 1], 1e-60)
   [1.10584379792120476018299547088585107443954623663875285836...]
   + [+/- ...]*I
   Connection between the regular singular points 0 and 1. The value at 1
   is the coordinate on the basis element of leading monomial 1.
   [ore_algebra] */
static void
test_fcc4_connection(void)
{
    gr_ctx_t ZZ, Pol, Dop, QQ;
    gr_ore_poly_t dop;
    gr_vec_t path;
    acb_t ref;
    acb_mat_t mat;
    slong prec = 256;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ctx_init_fmpq(QQ);
    gr_ore_poly_init(dop, Dop);
    gr_vec_init(path, 0, QQ);
    acb_init(ref);
    acb_mat_init(mat, 4, 4);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, FCC_DOP4, Dop));
    GR_MUST_SUCCEED(gr_vec_set_str(path, "[0, 1/2, 1]", 1, QQ));
    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path, QQ,
                                              NULL, NULL, 0, prec);
    check_status(status, "fcc4_connection");

    acb_set_str(ref,
        "[1.10584379792120476018299547088585107443954623663875285836"
        " +/- 1e-56] + [+/- 1e-48]*I", prec);
    check_acb(acb_mat_entry(mat, 0, 3), ref, prec - 20, "fcc4_connection");

    acb_mat_clear(mat);
    acb_clear(ref);
    gr_vec_clear(path, QQ);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(QQ);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* sage: dop.numerical_transition_matrix([75/17,4]).trace()
   [-2.048077478517672...e+95...] + [-2.818936813156118...e+95...]*I
   where dop is a random operator with prescribed singularities and exponents
   (from Chyzak, Goyer, Mezzarobba, 2022) that defines transition matrices
   with large entries, from the doctests of ore_algebra. Its leading
   coefficient is c*(x - 75/17)^3*(x - 4)^3. Both 75/17 and 4 are regular
   singular points. The default local bases at both ends are ordered in the
   reverse order compared to ore_algebra (by decreasing real part of the
   exponents), which leaves the trace unchanged. The linear solve currently
   loses about 200 bits.
   [ore_algebra] */
static void
test_large_entries(void)
{
    gr_ctx_t QQbar, Pol, Dop, QQ;
    gr_ore_poly_t dop;
    gr_vec_t path;
    acb_t ref, tr;
    acb_mat_t mat;
    slong prec = 384;
    int status = GR_SUCCESS;

    gr_ctx_init_complex_qqbar(QQbar);
    gr_ctx_init_gr_poly(Pol, QQbar);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ctx_init_fmpq(QQ);
    gr_ore_poly_init(dop, Dop);
    gr_vec_init(path, 0, QQ);
    acb_init(ref);
    acb_init(tr);
    acb_mat_init(mat, 3, 3);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop,
            "(21829048560*x^6 - 550862460720*x^5 + 5789380916880*x^4"
            " - 32434807101840*x^3 + 102165545592000*x^2"
            " - 171548863200000*x + 119964240000000)*D^3"
            " + (-506769209364*x^5 + 10798395588024*x^4"
            " - 91979496877620*x^3 + 391481933812992*x^2"
            " - 832564302926400*x + 707776470720000)*D^2"
            " + (251156588660*x^4 - 10005924091768*x^3"
            " + 99364591554260*x^2 - 379332696633792*x"
            " + 503580376843200)*D"
            " + 38987126730300*x^3 - 113173203800613*x^2"
            " - 1127681839519020*x + 3827710032458400", Dop));
    GR_MUST_SUCCEED(gr_vec_set_str(path, "[75/17, 143/34, 4]", 1, QQ));
    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path, QQ,
                                              NULL, NULL, 0, prec);
    check_status(status, "large_entries");
    acb_mat_trace(tr, mat, prec);

    acb_set_str(ref, "[-2.048077478517672e+95 +/- 1e+80]"
                     " + [-2.818936813156118e+95 +/- 1e+80]*I", prec);
    check_acb(tr, ref, prec - 216, "large_entries");

    acb_mat_clear(mat);
    acb_clear(tr);
    acb_clear(ref);
    gr_vec_clear(path, QQ);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(QQ);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(QQbar);
}

/* sage: dop = (x^2 - 2)*Dx^2 + x + 1
   sage: dop.numerical_transition_matrix([0, 1, QQbar(sqrt(2))], 1e-10)
   [         [2.49388146...] + [+/-...]*I          [2.40894178...] + [+/-...]*I]
   [[-0.203541775...] + [6.68738570...]*I  [0.204372067...] + [6.45961849...]*I]
   [ore_algebra] */
static void
test_sqrt2_connection(void)
{
    gr_ctx_t ZZ, Pol, Dop, QQbar, CC;
    gr_ore_poly_t dop;
    gr_vec_t path;
    acb_ptr ref;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_complex_qqbar(QQbar);
    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ctx_init_complex_acb(CC, prec);
    gr_ore_poly_init(dop, Dop);
    gr_vec_init(path, 0, QQbar);
    ref = _acb_vec_init(4);
    acb_mat_init(mat, 2, 2);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "(x^2 - 2)*D^2 + x + 1", Dop));
    GR_MUST_SUCCEED(gr_vec_set_str(path, "[0, 1, sqrt(2)]", 1, QQbar));
    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path, QQbar,
                                              NULL, NULL, 0, prec);
    check_status(status, "sqrt2_connection");

    GR_MUST_SUCCEED(_gr_vec_set_str(ref,
            "[[2.49388146 +/- 1e-8], [2.40894178 +/- 1e-8],"
            " [-0.203541775 +/- 1e-9] + [6.68738570 +/- 1e-8]*I,"
            " [0.204372067 +/- 1e-9] + [6.45961849 +/- 1e-8]*I]", 4, CC));
    for (slong i = 0; i < 2; i++)
        for (slong j = 0; j < 2; j++)
            check_acb(acb_mat_entry(mat, i, j), ref + 2 * i + j, prec - 12,
                      "sqrt2_connection");

    acb_mat_clear(mat);
    _acb_vec_clear(ref, 4);
    gr_vec_clear(path, QQbar);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(CC);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
    gr_ctx_clear(QQbar);
}

/* sage: _, z, Dz = DifferentialOperators(QQ, 'z')
   sage: dop = (z^2*(z^2-34*z+1)*Dz^3 + 3*z*(2*z^2-51*z+1)*Dz^2
   ....:       + (7*z^2-112*z+1)*Dz + (z-5))
   sage: roots = dop.leading_coefficient().roots(AA)
   sage: mat = dop.numerical_transition_matrix([0, roots[1][0]], 1e-10)
   sage: mat.list()
   [[4.846055616...], [-3.77845406...], [1.473024273...],
    [-14.9569783...]*I, [+/- ...]*I, [4.546376247...]*I,
    [-59.9006990...], [28.70759161...], [-18.2076291...]]
   sage: cst = -((1/4)*I)*(1+2^(1/2))^2*2^(3/4)/(pi*(2*2^(1/2)-3))
   sage: mat[1][2].overlaps(CBF(cst))
   True
   Both ends of the path are singular, and roots[1][0] = 17 - 12*sqrt(2).
   The local basis at that point is ordered differently from ore_algebra
   (rows 1 and 2 are swapped). The entry [+/- ...]*I is not checked.
   [ore_algebra] */
static void
test_apery_like(void)
{
    gr_ctx_t ZZ, Pol, Dop, QQbar, CC;
    gr_ore_poly_t dop;
    gr_vec_t path;
    acb_ptr ref;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_complex_qqbar(QQbar);
    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ctx_init_complex_acb(CC, prec);
    gr_ore_poly_init(dop, Dop);
    gr_vec_init(path, 0, QQbar);
    ref = _acb_vec_init(9);
    acb_mat_init(mat, 3, 3);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop,
            "x^2*(x^2 - 34*x + 1)*D^3 + 3*x*(2*x^2 - 51*x + 1)*D^2"
            " + (7*x^2 - 112*x + 1)*D + x - 5", Dop));
    GR_MUST_SUCCEED(gr_vec_set_str(path,
            "[0, 2/5*(17 - 12*sqrt(2)), 3/5*(17 - 12*sqrt(2)),"
            " 17 - 12*sqrt(2)]", 1, QQbar));
    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path, QQbar,
                                              NULL, NULL, 0, prec);
    check_status(status, "apery_like");

    GR_MUST_SUCCEED(_gr_vec_set_str(ref,
            "[[4.846055616 +/- 1e-9], [-3.77845406 +/- 1e-8],"
            " [1.473024273 +/- 1e-9],"
            " [-59.9006990 +/- 1e-7], [28.70759161 +/- 1e-8],"
            " [-18.2076291 +/- 1e-7],"
            " [-14.9569783 +/- 1e-7]*I, 0,"
            " -I/4*(1 + sqrt(2))^2*2^(3/4)/(pi*(2*sqrt(2) - 3))]", 9, CC));
    for (slong i = 0; i < 3; i++)
        for (slong j = 0; j < 3; j++)
            if (i != 2 || j != 1)
                check_acb(acb_mat_entry(mat, i, j), ref + 3 * i + j,
                          prec - 16, "apery_like");

    acb_mat_clear(mat);
    _acb_vec_clear(ref, 9);
    gr_vec_clear(path, QQbar);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(CC);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
    gr_ctx_clear(QQbar);
}

/* sage: dop = (x^2 + 1)*Dx^2 + 2*x*Dx
   sage: dop.numerical_transition_matrix([0, i + 1, 2*i, i - 1, 0])
   [ [1.00...] + [+/- ...]*I  [3.141592653589793...] + [+/-...]*I]
   [ [+/- ...] + [+/- ...]*I  [1.000000000000000...] + [+/-...]*I]
   Monodromy of arctan around i, along an octagon of radius 1/2.
   [ore_algebra] */
static void
test_arctan_monodromy(void)
{
    gr_ctx_t ZZ, Pol, Dop, CC;
    gr_ore_poly_t dop;
    gr_vec_t path;
    acb_ptr ref;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ctx_init_complex_acb(CC, prec);
    gr_ore_poly_init(dop, Dop);
    gr_vec_init(path, 0, CC);
    ref = _acb_vec_init(4);
    acb_mat_init(mat, 2, 2);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "(x^2 + 1)*D^2 + 2*x*D", Dop));
    GR_MUST_SUCCEED(gr_vec_set_str(path,
            "[0, I/2, 3/8 + 5/8*I, 1/2 + I, 3/8 + 11/8*I, 3/2*I,"
            " -3/8 + 11/8*I, -1/2 + I, -3/8 + 5/8*I, I/2, I/4, 0]", 1, CC));
    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path, CC,
                                              NULL, NULL, 0, prec);
    check_status(status, "arctan_monodromy");

    GR_MUST_SUCCEED(_gr_vec_set_str(ref, "[1, pi, 0, 1]", 4, CC));
    for (slong i = 0; i < 2; i++)
        for (slong j = 0; j < 2; j++)
            check_acb(acb_mat_entry(mat, i, j), ref + 2 * i + j, prec - 12,
                      "arctan_monodromy");

    acb_mat_clear(mat);
    _acb_vec_clear(ref, 4);
    gr_vec_clear(path, CC);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(CC);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
}

/* sage: dop.numerical_transition_matrix([0, i], 1e-10)
   [            ...           [+/-...] + [-0.50000000...]*I]
   [ ...1.000000...  [0.7853981634...] + [0.346573590...]*I]
   sage: dop.numerical_transition_matrix([0, i - 1, i], 1e-10)
   [     [+/-...] + [+/-...]*I         [+/-...] + [-0.5000000000...]*I]
   [ [1.00000...] + [+/-...]*I [-2.356194490...] + [0.3465735902...]*I]
   (with dop as in test_arctan_monodromy) Near i, one has
   arctan(x) = -i/2*log(x - i) + pi/4 + i*log(2)/2 + O(x - i), and the
   second path ends on the branch cut, coming from below.
   [ore_algebra] */
static void
test_arctan_connection(void)
{
    gr_ctx_t ZZ, Pol, Dop, QQbar, CC;
    gr_ore_poly_t dop;
    gr_vec_t path;
    acb_ptr ref;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_complex_qqbar(QQbar);
    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ctx_init_complex_acb(CC, prec);
    gr_ore_poly_init(dop, Dop);
    gr_vec_init(path, 0, QQbar);
    ref = _acb_vec_init(4);
    acb_mat_init(mat, 2, 2);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "(x^2 + 1)*D^2 + 2*x*D", Dop));

    GR_MUST_SUCCEED(gr_vec_set_str(path, "[0, I]", 1, QQbar));
    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path, QQbar,
                                              NULL, NULL, 0, prec);
    check_status(status, "arctan_connection");

    GR_MUST_SUCCEED(_gr_vec_set_str(ref,
            "[0, -I/2, 1, pi/4 + I*log(2)/2]", 4, CC));
    for (slong i = 0; i < 2; i++)
        for (slong j = 0; j < 2; j++)
            check_acb(acb_mat_entry(mat, i, j), ref + 2 * i + j, prec - 8,
                      "arctan_connection");

    GR_MUST_SUCCEED(gr_vec_set_str(path, "[0, -3/4 + I/4, -1 + I, I]", 1,
                                   QQbar));
    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path, QQbar,
                                              NULL, NULL, 0, prec);
    check_status(status, "arctan_connection");

    GR_MUST_SUCCEED(_gr_vec_set_str(ref,
            "[0, -I/2, 1, -3*pi/4 + I*log(2)/2]", 4, CC));
    for (slong i = 0; i < 2; i++)
        for (slong j = 0; j < 2; j++)
            check_acb(acb_mat_entry(mat, i, j), ref + 2 * i + j, prec - 12,
                      "arctan_connection");

    acb_mat_clear(mat);
    _acb_vec_clear(ref, 4);
    gr_vec_clear(path, QQbar);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(CC);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
    gr_ctx_clear(QQbar);
}

/* sage: dop = x*Dx^3 + 2*Dx^2 + x*Dx
   sage: mat = dop.numerical_transition_matrix([-1, 0, i, -1])
   sage: id = identity_matrix(3)
   sage: all(y.rad() < 1e-13 for row in (mat - id) for y in row)
   True
   The path passes through the singular point 0, arriving along the branch
   cut and leaving upwards, so that it does not wind around 0.
   [ore_algebra] */
static void
test_loop_through_singularity(void)
{
    gr_ctx_t ZZ, Pol, Dop, QQbar;
    gr_ore_poly_t dop;
    gr_vec_t path;
    acb_t c;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_complex_qqbar(QQbar);
    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    gr_vec_init(path, 0, QQbar);
    acb_init(c);
    acb_mat_init(mat, 3, 3);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "x*D^3 + 2*D^2 + x*D", Dop));
    GR_MUST_SUCCEED(gr_vec_set_str(path, "[-1, 0, I, -3/4 + 3/4*I, -1]", 1,
                                   QQbar));
    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path, QQbar,
                                              NULL, NULL, 0, prec);
    check_status(status, "loop_through_singularity");

    for (slong i = 0; i < 3; i++)
    {
        for (slong j = 0; j < 3; j++)
        {
            acb_set_si(c, i == j);
            check_acb(acb_mat_entry(mat, i, j), c, prec - 12,
                      "loop_through_singularity");
        }
    }

    acb_mat_clear(mat);
    acb_clear(c);
    gr_vec_clear(path, QQbar);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
    gr_ctx_clear(QQbar);
}

/* sage: a = QQbar(2*exp(I*2*pi/3))
   sage: b = a + a.conjugate()
   sage: ((x-1)*Dx*x*Dx).numerical_transition_matrix([0, b])[0,0]
   [0.6931471805599453...] + [3.141592653589793...]*I
   sage: ((x-1)*Dx*x*Dx).numerical_transition_matrix([0, b - I/10^1000])[0,0]
   [0.6931471805599453...] + [-3.141592653589793...]*I
   Here b = -2, and the entry is the value at b of the solution log(x). The
   second path leaves 0 through the lower half-plane instead of perturbing b
   (which would make the exact Taylor shifts very expensive). The step from
   -1/2 - i/4 to -9/10 currently loses about 10% of the working precision.
   [ore_algebra] */
static void
test_log_subdivided(void)
{
    gr_ctx_t ZZ, Pol, Dop, QQbar, QQ;
    gr_ore_poly_t dop;
    gr_vec_t path0, path1;
    acb_t ref;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_complex_qqbar(QQbar);
    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ctx_init_fmpq(QQ);
    gr_ore_poly_init(dop, Dop);
    gr_vec_init(path0, 0, QQ);
    gr_vec_init(path1, 0, QQbar);
    acb_init(ref);
    acb_mat_init(mat, 2, 2);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "(x^2 - x)*D^2 + (x - 1)*D",
                                        Dop));

    GR_MUST_SUCCEED(gr_vec_set_str(path0, "[0, -1/2, -9/10, -8/5, -2]", 1,
                                   QQ));
    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path0, QQ,
                                              NULL, NULL, 0, prec);
    check_status(status, "log_subdivided");
    acb_set_str(ref, "log(2) + pi*I", prec);
    check_acb(acb_mat_entry(mat, 0, 0), ref, prec - 4, "log_subdivided");

    GR_MUST_SUCCEED(gr_vec_set_str(path1,
            "[0, -1/2 - I/4, -9/10, -8/5, -2]", 1, QQbar));
    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path1, QQbar,
                                              NULL, NULL, 0, prec);
    check_status(status, "log_subdivided");
    acb_set_str(ref, "log(2) - pi*I", prec);
    check_acb(acb_mat_entry(mat, 0, 0), ref, prec - 24, "log_subdivided");

    acb_mat_clear(mat);
    acb_clear(ref);
    gr_vec_clear(path1, QQbar);
    gr_vec_clear(path0, QQ);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(QQ);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);
    gr_ctx_clear(QQbar);
}

/* sage: ((9*x^2 - 1)*Dx + 18*x).numerical_transition_matrix([0,I,1],
   ....:         algorithm="naive", squash_intervals=True)
   [[-0.125000000000000...] + [+/- ...]*I]
   The solutions are the multiples of 1/(9*x^2 - 1).
   [ore_algebra] */
static void
test_rational_solution(void)
{
    gr_ctx_t CC, Pol, Dop;
    gr_ore_poly_t dop;
    gr_vec_t path;
    acb_t ref;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_complex_acb(CC, prec);
    gr_ctx_init_gr_poly(Pol, CC);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ore_poly_init(dop, Dop);
    gr_vec_init(path, 0, CC);
    acb_init(ref);
    acb_mat_init(mat, 1, 1);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop, "(9*x^2 - 1)*D + 18*x", Dop));
    GR_MUST_SUCCEED(gr_vec_set_str(path,
            "[0, I/4, I/2, I, 1/2 + 3/4*I, 1 + I/2, 1]", 1, CC));
    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path, CC,
                                              NULL, NULL, 0, prec);
    check_status(status, "rational_solution");

    acb_set_str(ref, "-1/8", prec);
    check_acb(acb_mat_entry(mat, 0, 0), ref, prec - 8, "rational_solution");

    acb_mat_clear(mat);
    acb_clear(ref);
    gr_vec_clear(path, CC);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(CC);
}

/* sage: dop =  -452*Dx^10 + (-2*x^2 - x - 1/2)*Dx^9 + ...
   sage: ((dop.numerical_transition_matrix([0,1])*dop.numerical_transition_matrix([1, 1+i, 0]) - 1)).norm('frob') < 1e-13
   True
   The operator has no singular points.
   [ore_algebra] */
static void
test_entire_loop(void)
{
    gr_ctx_t QQ, Pol, Dop, CC;
    gr_ore_poly_t dop;
    gr_vec_t path;
    acb_t c;
    acb_mat_t mat;
    slong prec = 128;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_gr_poly(Pol, QQ);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    gr_ctx_init_complex_acb(CC, prec);
    gr_ore_poly_init(dop, Dop);
    gr_vec_init(path, 0, CC);
    acb_init(c);
    acb_mat_init(mat, 10, 10);

    GR_MUST_SUCCEED(gr_ore_poly_set_str(dop,
            "-452*D^10 + (-2*x^2 - x - 1/2)*D^9 + (1/2*x + 22)*D^8"
            " + (1/4*x^2 + x)*D^7 + (1/3*x^2 - 1/2*x + 1/3)*D^6"
            " + (-3*x^2 + x + 1)*D^5 + (x^2 - 4/3*x)*D^4 + (2*x^2 - 2*x)*D^3"
            " + (2*x^2 + x)*D^2 + (-2/3*x^2 - 5/27*x - 1/3)*D"
            " - 18*x^2 + 6/5*x - 6", Dop));
    GR_MUST_SUCCEED(gr_vec_set_str(path, "[0, 1, 1 + I, 0]", 1, CC));
    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path, CC,
                                              NULL, NULL, 0, prec);
    check_status(status, "entire_loop");

    for (slong i = 0; i < 10; i++)
    {
        for (slong j = 0; j < 10; j++)
        {
            acb_set_si(c, i == j);
            check_acb(acb_mat_entry(mat, i, j), c, prec - 24, "entire_loop");
        }
    }

    acb_mat_clear(mat);
    acb_clear(c);
    gr_vec_clear(path, CC);
    gr_ore_poly_clear(dop, Dop);
    gr_ctx_clear(CC);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(QQ);
}

TEST_FUNCTION_START(acb_ode_transition_matrix_fixed_vs_ref, state)
{
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
    test_sqrt2_connection();
    test_apery_like();
    test_arctan_monodromy();
    test_arctan_connection();
    test_loop_through_singularity();
    test_log_subdivided();
    test_rational_solution();
    test_entire_loop();

    TEST_FUNCTION_END(state);
}
