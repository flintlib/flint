/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    A catalog of test inputs for the Bessel functions, the incomplete gamma
    functions and the hypergeometric functions 1F1, U and 2F1, with the
    time and the accuracy of the default functions and, optionally, of the
    separate algorithms and of loops which increase the precision.

    Build with "make profile MOD=acb_hypgeom" and run

        build/acb_hypgeom/profile/p-catalog -help

    for the options.  Some examples:

        p-catalog                               everything at 64, 333, 3333 bits
        p-catalog -list                         the catalog with its classes
        p-catalog -func bessel -prec 64         J, Y, I, K at 64 bits
        p-catalog -func 2f1 -class plarge+osc   one class of inputs
        p-catalog -func bessel_j -alg all       all algorithms for J
        p-catalog -alg all,loop,loop:nointegration -ref -summary -csv out.csv
        p-catalog -case K41 -prec 333 -alg all -ref -values
        p-catalog -input "2f1(5, 2, 400, 11/10*i)" -alg all,loop

    The inputs are expressions for gr_set_str which are evaluated at the
    working precision; "prec" in an expression stands for that precision,
    and "x +/- r" for a ball.  The catalog is the set of tables below
    ("The catalog").  To keep the identifiers of the cases, new cases are
    appended at the end of a table.

    Each measurement is made in a child process where fork is available,
    so that an algorithm which is forced on unsuitable input can be
    stopped when it exceeds the time limit, the memory limit, or aborts.
*/

/* fork, poll, kill, setrlimit */
#if !defined(_WIN32) && !defined(_MSC_VER) && !defined(_POSIX_C_SOURCE)
# define _POSIX_C_SOURCE 200809L
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "profiler.h"
#include "arb.h"
#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "acb_hypgeom/impl.h"
#include "gr.h"

#if !defined(CATALOG_NO_FORK) \
    && ((defined(__unix__) && !defined(__CYGWIN__) && !defined(__EMSCRIPTEN__)) || defined(__APPLE__))
# define CATALOG_HAVE_FORK 1
# include <errno.h>
# include <fcntl.h>
# include <poll.h>
# include <signal.h>
# include <unistd.h>
# include <sys/resource.h>
# include <sys/wait.h>
#else
# define CATALOG_HAVE_FORK 0
#endif

#define MAX_ARGS 4
#define MAX_PRECS 32
#define MAX_ALGS 16
#define MAX_LOOP_STEPS 24
#define MAX_REF_POINTS (1 + 4 * MAX_ARGS)
#define MAX_INPUTS 64

/* ------------------------------------------------------------------------- */
/*  Functions and algorithms                                                 */
/* ------------------------------------------------------------------------- */

/* An algorithm evaluates the function at x[0], ..., x[nargs-1]; it returns
   0 if it does not apply to this input (and 1 otherwise, with a result
   which is not finite if the algorithm fails). */
typedef int (*eval_func)(acb_t res, acb_srcptr x, slong prec);

typedef struct
{
    const char * name;
    eval_func f;
    const char * desc;
}
alg_struct;

typedef struct
{
    const char * name;
    char id;                    /* the letter of the case identifiers */
    int nargs;
    const char * argnames[MAX_ARGS];
    const char * desc;
    const alg_struct * algs;    /* the first is the default function */
}
func_struct;

#define ALG(name, body) \
    static int name(acb_t res, acb_srcptr x, slong prec) \
    { body; return 1; }

/* The functions of arb_hypgeom, for real input */
#define ALG_REAL(name, n, body) \
    static int name(acb_t res, acb_srcptr x, slong prec) \
    { \
        slong i; \
        for (i = 0; i < n; i++) \
            if (!acb_is_real(x + i)) \
                return 0; \
        body; \
        arb_zero(acb_imagref(res)); \
        return 1; \
    }

#define RES acb_realref(res)
#define RE(i) acb_realref(x + (i))

/* Bessel functions */

ALG(j_default, acb_hypgeom_bessel_j(res, x, x + 1, prec))
ALG(j_0f1, acb_hypgeom_bessel_j_0f1(res, x, x + 1, prec))
ALG(j_asymp, acb_hypgeom_bessel_j_asymp(res, x, x + 1, prec))
ALG(j_nointegration, acb_hypgeom_bessel_j_nointegration(res, x, x + 1, prec))
ALG(j_integration, acb_hypgeom_bessel_j_integration(res, x, x + 1, prec))
ALG_REAL(j_arb, 2, arb_hypgeom_bessel_j(RES, RE(0), RE(1), prec))

ALG(y_default, acb_hypgeom_bessel_y(res, x, x + 1, prec))
ALG(y_nointegration, acb_hypgeom_bessel_y_nointegration(res, x, x + 1, prec))
ALG(y_integration, acb_hypgeom_bessel_y_integration(res, x, x + 1, prec))
ALG_REAL(y_arb, 2, arb_hypgeom_bessel_y(RES, RE(0), RE(1), prec))

ALG(i_default, acb_hypgeom_bessel_i(res, x, x + 1, prec))
ALG(i_0f1, acb_hypgeom_bessel_i_0f1(res, x, x + 1, 0, prec))
ALG(i_asymp, acb_hypgeom_bessel_i_asymp(res, x, x + 1, 0, prec))
ALG(i_nointegration, acb_hypgeom_bessel_i_nointegration(res, x, x + 1, prec))
ALG(i_integration, acb_hypgeom_bessel_i_integration(res, x, x + 1, prec))
ALG_REAL(i_arb, 2, arb_hypgeom_bessel_i(RES, RE(0), RE(1), prec))

ALG(k_default, acb_hypgeom_bessel_k(res, x, x + 1, prec))
ALG(k_0f1, acb_hypgeom_bessel_k_0f1(res, x, x + 1, 0, prec))
ALG(k_asymp, acb_hypgeom_bessel_k_asymp(res, x, x + 1, 0, prec))
ALG(k_nointegration, acb_hypgeom_bessel_k_nointegration(res, x, x + 1, prec))
ALG(k_integration, acb_hypgeom_bessel_k_integration(res, x, x + 1, prec))
ALG_REAL(k_arb, 2, arb_hypgeom_bessel_k(RES, RE(0), RE(1), prec))

static const alg_struct algs_j[] = {
    { "default", j_default, "acb_hypgeom_bessel_j" },
    { "0f1", j_0f1, "acb_hypgeom_bessel_j_0f1 (convergent series)" },
    { "asymp", j_asymp, "acb_hypgeom_bessel_j_asymp (asymptotic expansion)" },
    { "nointegration", j_nointegration, "acb_hypgeom_bessel_j_nointegration (series or asymptotic expansion)" },
    { "integration", j_integration, "acb_hypgeom_bessel_j_integration (steepest descent)" },
    { "arb", j_arb, "arb_hypgeom_bessel_j (real input)" },
    { NULL, NULL, NULL }
};

static const alg_struct algs_y[] = {
    { "default", y_default, "acb_hypgeom_bessel_y" },
    { "nointegration", y_nointegration, "acb_hypgeom_bessel_y_nointegration (series or asymptotic expansion)" },
    { "integration", y_integration, "acb_hypgeom_bessel_y_integration (steepest descent)" },
    { "arb", y_arb, "arb_hypgeom_bessel_y (real input)" },
    { NULL, NULL, NULL }
};

static const alg_struct algs_i[] = {
    { "default", i_default, "acb_hypgeom_bessel_i" },
    { "0f1", i_0f1, "acb_hypgeom_bessel_i_0f1 (convergent series)" },
    { "asymp", i_asymp, "acb_hypgeom_bessel_i_asymp (asymptotic expansion)" },
    { "nointegration", i_nointegration, "acb_hypgeom_bessel_i_nointegration (series or asymptotic expansion)" },
    { "integration", i_integration, "acb_hypgeom_bessel_i_integration (steepest descent)" },
    { "arb", i_arb, "arb_hypgeom_bessel_i (real input)" },
    { NULL, NULL, NULL }
};

static const alg_struct algs_k[] = {
    { "default", k_default, "acb_hypgeom_bessel_k" },
    { "0f1", k_0f1, "acb_hypgeom_bessel_k_0f1 (convergent series)" },
    { "asymp", k_asymp, "acb_hypgeom_bessel_k_asymp (asymptotic expansion)" },
    { "nointegration", k_nointegration, "acb_hypgeom_bessel_k_nointegration (series or asymptotic expansion)" },
    { "integration", k_integration, "acb_hypgeom_bessel_k_integration (steepest descent)" },
    { "arb", k_arb, "arb_hypgeom_bessel_k (real input)" },
    { NULL, NULL, NULL }
};

/* Incomplete gamma functions */

static int
_gamma_integration(acb_t res, acb_srcptr x, int upper, slong prec)
{
    if (!_acb_hypgeom_gamma_incomplete_integration(res, x, x + 1, upper, prec))
        acb_indeterminate(res);
    return 1;
}

static int
_gamma_singular(acb_t res, acb_srcptr x, int regularized, slong prec)
{
    if (!acb_is_int(x) || arf_sgn(arb_midref(acb_realref(x))) > 0
        || arf_cmpabs_2exp_si(arb_midref(acb_realref(x)), 30) >= 0)
        return 0;

    acb_hypgeom_gamma_upper_singular(res,
        arf_get_si(arb_midref(acb_realref(x)), ARF_RND_DOWN), x + 1, regularized, prec);
    return 1;
}

ALG(gu_default, acb_hypgeom_gamma_upper(res, x, x + 1, 0, prec))
ALG(gu_asymp, acb_hypgeom_gamma_upper_asymp(res, x, x + 1, 0, prec))
ALG(gu_1f1a, acb_hypgeom_gamma_upper_1f1a(res, x, x + 1, 0, prec))
ALG(gu_1f1b, acb_hypgeom_gamma_upper_1f1b(res, x, x + 1, 0, prec))
ALG(gu_nointegration, acb_hypgeom_gamma_upper_nointegration(res, x, x + 1, 0, prec))
ALG_REAL(gu_arb, 2, arb_hypgeom_gamma_upper(RES, RE(0), RE(1), 0, prec))
static int gu_singular(acb_t res, acb_srcptr x, slong prec) { return _gamma_singular(res, x, 0, prec); }
static int gu_integration(acb_t res, acb_srcptr x, slong prec) { return _gamma_integration(res, x, 1, prec); }

ALG(gl_default, acb_hypgeom_gamma_lower(res, x, x + 1, 0, prec))
ALG(gl_nointegration, acb_hypgeom_gamma_lower_nointegration(res, x, x + 1, 0, prec))
ALG_REAL(gl_arb, 2, arb_hypgeom_gamma_lower(RES, RE(0), RE(1), 0, prec))
static int gl_integration(acb_t res, acb_srcptr x, slong prec) { return _gamma_integration(res, x, 0, prec); }

ALG(gq_default, acb_hypgeom_gamma_upper(res, x, x + 1, 1, prec))
ALG(gq_asymp, acb_hypgeom_gamma_upper_asymp(res, x, x + 1, 1, prec))
ALG(gq_1f1a, acb_hypgeom_gamma_upper_1f1a(res, x, x + 1, 1, prec))
ALG(gq_1f1b, acb_hypgeom_gamma_upper_1f1b(res, x, x + 1, 1, prec))
ALG(gq_nointegration, acb_hypgeom_gamma_upper_nointegration(res, x, x + 1, 1, prec))
ALG_REAL(gq_arb, 2, arb_hypgeom_gamma_upper(RES, RE(0), RE(1), 1, prec))

ALG(gp_default, acb_hypgeom_gamma_lower(res, x, x + 1, 1, prec))
ALG(gp_nointegration, acb_hypgeom_gamma_lower_nointegration(res, x, x + 1, 1, prec))
ALG_REAL(gp_arb, 2, arb_hypgeom_gamma_lower(RES, RE(0), RE(1), 1, prec))

static const alg_struct algs_gu[] = {
    { "default", gu_default, "acb_hypgeom_gamma_upper" },
    { "asymp", gu_asymp, "acb_hypgeom_gamma_upper_asymp (asymptotic expansion)" },
    { "1f1a", gu_1f1a, "acb_hypgeom_gamma_upper_1f1a (series for the lower function)" },
    { "1f1b", gu_1f1b, "acb_hypgeom_gamma_upper_1f1b (series, Kummer's transformation)" },
    { "singular", gu_singular, "acb_hypgeom_gamma_upper_singular (s = 0, -1, -2, ...)" },
    { "nointegration", gu_nointegration, "acb_hypgeom_gamma_upper_nointegration" },
    { "integration", gu_integration, "_acb_hypgeom_gamma_incomplete_integration (steepest descent)" },
    { "arb", gu_arb, "arb_hypgeom_gamma_upper (real input)" },
    { NULL, NULL, NULL }
};

static const alg_struct algs_gl[] = {
    { "default", gl_default, "acb_hypgeom_gamma_lower" },
    { "nointegration", gl_nointegration, "acb_hypgeom_gamma_lower_nointegration" },
    { "integration", gl_integration, "_acb_hypgeom_gamma_incomplete_integration (steepest descent)" },
    { "arb", gl_arb, "arb_hypgeom_gamma_lower (real input)" },
    { NULL, NULL, NULL }
};

static const alg_struct algs_gq[] = {
    { "default", gq_default, "acb_hypgeom_gamma_upper, regularized" },
    { "asymp", gq_asymp, "acb_hypgeom_gamma_upper_asymp (asymptotic expansion)" },
    { "1f1a", gq_1f1a, "acb_hypgeom_gamma_upper_1f1a (series for the lower function)" },
    { "1f1b", gq_1f1b, "acb_hypgeom_gamma_upper_1f1b (series, Kummer's transformation)" },
    { "nointegration", gq_nointegration, "acb_hypgeom_gamma_upper_nointegration" },
    { "arb", gq_arb, "arb_hypgeom_gamma_upper (real input)" },
    { NULL, NULL, NULL }
};

static const alg_struct algs_gp[] = {
    { "default", gp_default, "acb_hypgeom_gamma_lower, regularized" },
    { "nointegration", gp_nointegration, "acb_hypgeom_gamma_lower_nointegration" },
    { "arb", gp_arb, "arb_hypgeom_gamma_lower (real input)" },
    { NULL, NULL, NULL }
};

/* 1F1 and U */

static int
u_asymp(acb_t res, acb_srcptr x, slong prec)
{
    acb_t t;
    acb_init(t);
    /* U = z^(-a) U* */
    acb_hypgeom_u_asymp(res, x, x + 1, x + 2, -1, prec);
    acb_neg(t, x);
    acb_pow(t, x + 2, t, prec);
    acb_mul(res, res, t, prec);
    acb_clear(t);
    return 1;
}

ALG(m_default, acb_hypgeom_m(res, x, x + 1, x + 2, 0, prec))
ALG(m_1f1, acb_hypgeom_m_1f1(res, x, x + 1, x + 2, 0, prec))
ALG(m_asymp, acb_hypgeom_m_asymp(res, x, x + 1, x + 2, 0, prec))
ALG(m_nointegration, _acb_hypgeom_m_nointegration(res, x, x + 1, x + 2, 0, prec))
ALG(m_integration, acb_hypgeom_1f1_integration(res, x, x + 1, x + 2, 0, prec))
ALG_REAL(m_arb, 3, arb_hypgeom_m(RES, RE(0), RE(1), RE(2), 0, prec))

ALG(u_default, acb_hypgeom_u(res, x, x + 1, x + 2, prec))
ALG(u_1f1, acb_hypgeom_u_1f1(res, x, x + 1, x + 2, prec))
ALG(u_nointegration, _acb_hypgeom_u_nointegration(res, x, x + 1, x + 2, prec))
ALG(u_integration, acb_hypgeom_u_integration(res, x, x + 1, x + 2, prec))
ALG_REAL(u_arb, 3, arb_hypgeom_u(RES, RE(0), RE(1), RE(2), prec))

static const alg_struct algs_m[] = {
    { "default", m_default, "acb_hypgeom_m" },
    { "1f1", m_1f1, "acb_hypgeom_m_1f1 (convergent series)" },
    { "asymp", m_asymp, "acb_hypgeom_m_asymp (asymptotic expansion)" },
    { "nointegration", m_nointegration, "_acb_hypgeom_m_nointegration (series or asymptotic expansion)" },
    { "integration", m_integration, "acb_hypgeom_1f1_integration (steepest descent)" },
    { "arb", m_arb, "arb_hypgeom_m (real input)" },
    { NULL, NULL, NULL }
};

static const alg_struct algs_u[] = {
    { "default", u_default, "acb_hypgeom_u" },
    { "1f1", u_1f1, "acb_hypgeom_u_1f1 (convergent series)" },
    { "asymp", u_asymp, "acb_hypgeom_u_asymp (asymptotic expansion)" },
    { "nointegration", u_nointegration, "_acb_hypgeom_u_nointegration (series or asymptotic expansion)" },
    { "integration", u_integration, "acb_hypgeom_u_integration (steepest descent)" },
    { "arb", u_arb, "arb_hypgeom_u (real input)" },
    { NULL, NULL, NULL }
};

/* 2F1 */

ALG(f_default, acb_hypgeom_2f1(res, x, x + 1, x + 2, x + 3, 0, prec))
ALG(f_direct, acb_hypgeom_2f1_direct(res, x, x + 1, x + 2, x + 3, 0, prec))
ALG(f_transform1, acb_hypgeom_2f1_transform(res, x, x + 1, x + 2, x + 3, 0, 1, prec))
ALG(f_transform2, acb_hypgeom_2f1_transform(res, x, x + 1, x + 2, x + 3, 0, 2, prec))
ALG(f_transform3, acb_hypgeom_2f1_transform(res, x, x + 1, x + 2, x + 3, 0, 3, prec))
ALG(f_transform4, acb_hypgeom_2f1_transform(res, x, x + 1, x + 2, x + 3, 0, 4, prec))
ALG(f_transform5, acb_hypgeom_2f1_transform(res, x, x + 1, x + 2, x + 3, 0, 5, prec))
ALG(f_corner, acb_hypgeom_2f1_corner(res, x, x + 1, x + 2, x + 3, 0, prec))
ALG(f_nointegration, _acb_hypgeom_2f1_nointegration(res, x, x + 1, x + 2, x + 3, 0, prec))
ALG(f_integration, acb_hypgeom_2f1_integration(res, x, x + 1, x + 2, x + 3, 0, prec))
ALG_REAL(f_arb, 4, arb_hypgeom_2f1(RES, RE(0), RE(1), RE(2), RE(3), 0, prec))

static const alg_struct algs_f[] = {
    { "default", f_default, "acb_hypgeom_2f1" },
    { "direct", f_direct, "acb_hypgeom_2f1_direct (series at z)" },
    { "transform1", f_transform1, "acb_hypgeom_2f1_transform, z/(z-1)" },
    { "transform2", f_transform2, "acb_hypgeom_2f1_transform, 1/z" },
    { "transform3", f_transform3, "acb_hypgeom_2f1_transform, 1/(1-z)" },
    { "transform4", f_transform4, "acb_hypgeom_2f1_transform, 1-z" },
    { "transform5", f_transform5, "acb_hypgeom_2f1_transform, 1-1/z" },
    { "corner", f_corner, "acb_hypgeom_2f1_corner (z near exp(+-pi i/3))" },
    { "nointegration", f_nointegration, "_acb_hypgeom_2f1_nointegration (method chosen from z)" },
    { "integration", f_integration, "acb_hypgeom_2f1_integration (steepest descent)" },
    { "arb", f_arb, "arb_hypgeom_2f1 (real input)" },
    { NULL, NULL, NULL }
};

static const func_struct funcs[] = {
    { "bessel_j", 'J', 2, { "nu", "z" }, "Bessel function J_nu(z)", algs_j },
    { "bessel_y", 'Y', 2, { "nu", "z" }, "Bessel function Y_nu(z)", algs_y },
    { "bessel_i", 'I', 2, { "nu", "z" }, "modified Bessel function I_nu(z)", algs_i },
    { "bessel_k", 'K', 2, { "nu", "z" }, "modified Bessel function K_nu(z)", algs_k },
    { "gamma_upper", 'G', 2, { "s", "z" }, "upper incomplete gamma function Gamma(s,z)", algs_gu },
    { "gamma_lower", 'L', 2, { "s", "z" }, "lower incomplete gamma function gamma(s,z)", algs_gl },
    { "gamma_upper_reg", 'Q', 2, { "s", "z" }, "Q(s,z) = Gamma(s,z) / Gamma(s)", algs_gq },
    { "gamma_lower_reg", 'P', 2, { "s", "z" }, "P(s,z) = gamma(s,z) / Gamma(s)", algs_gp },
    { "1f1", 'M', 3, { "a", "b", "z" }, "confluent hypergeometric function 1F1(a,b,z) = M(a,b,z)", algs_m },
    { "u", 'U', 3, { "a", "b", "z" }, "confluent hypergeometric function U(a,b,z)", algs_u },
    { "2f1", 'F', 4, { "a", "b", "c", "z" }, "hypergeometric function 2F1(a,b,c,z)", algs_f },
};

#define NUM_FUNCS ((slong) (sizeof(funcs) / sizeof(funcs[0])))

/* ------------------------------------------------------------------------- */
/*  The catalog                                                              */
/* ------------------------------------------------------------------------- */

/*
    An entry gives the functions to which it applies (their letters, or
    NULL for the default of the table), the inputs, and the classes which
    are not found automatically:

    osc         oscillatory region
    mono        exponentially small or large (monotonic) region
    transition  argument comparable to the parameters
    turning     at or near a turning point
    asymp       argument much larger than the parameters
    series      argument much smaller than the parameters
    zero        near a zero of the function
    cancel      cancellation (inherent, or in the series)
    cut         on or near a branch cut
    corner      2F1 near exp(+-pi i/3)
    poly        terminating series
    degenerate  at or near a limit case of a connection formula
    pole        at or near a pole in a parameter
    special     reduces to a simpler function
    wide        wide balls

    The classes found from the inputs are exact, rounded (inexact from the
    evaluation of the expression only), inexact (with a radius, which is
    fixedrad or precrad), real, complex, int, halfint, nearint (some
    parameter), ptiny, pnormal, plarge, phuge (the largest parameter;
    ptiny also for any tiny nonzero parameter) and ztiny, znormal, zlarge,
    zhuge (the argument), where tiny < 2^-10 <= normal < 64 <= large
    < 30000 <= huge.
*/

typedef struct
{
    const char * sel;
    const char * args[MAX_ARGS];
    const char * tags;
}
entry_struct;

#define PR "(1 +/- 2^-(prec-8))"    /* accurate to 8 bits less than the precision */

/* J_nu(z) and Y_nu(z) */
static const entry_struct entries_jy[] = {
    /* orders and arguments of moderate size */
    { NULL, { "0", "1" }, "" },
    { NULL, { "1", "1/2" }, "" },
    { NULL, { "0", "10" }, "osc" },
    { NULL, { "1/2", "3" }, "special" },
    { NULL, { "7/2", "10" }, "special" },
    { NULL, { "1/3", "5" }, "" },
    { NULL, { "-1/3", "5" }, "" },
    { NULL, { "-5/2", "2" }, "special" },
    { NULL, { "-3", "2" }, "" },
    { NULL, { "5", "5" }, "" },
    { NULL, { "pi", "sqrt(2)" }, "" },
    /* tiny arguments and orders */
    { NULL, { "0", "10^-30" }, "series" },
    { NULL, { "1", "10^-30" }, "series" },
    { NULL, { "1/3", "10^-100" }, "series" },
    { NULL, { "10^-20", "1" }, "degenerate" },
    { NULL, { "10^-30", "10^-30" }, "degenerate series" },
    /* orders near integers and half-integers */
    { NULL, { "2 + 2^-30", "3" }, "degenerate" },
    { NULL, { "2 + 10^-15", "3" }, "degenerate" },
    { NULL, { "-3 + 2^-40", "2" }, "degenerate" },
    { NULL, { "1/2 + 2^-40", "3" }, "" },
    /* near zeros */
    { "J", { "0", "2.4048255576957727686" }, "zero cancel" },
    { "J", { "10", "14.475500686554541238451637655413" }, "zero cancel" },
    { "Y", { "0", "0.89357696627916752158" }, "zero cancel" },
    { "Y", { "10", "12.128927704415439386692455063522" }, "zero cancel" },
    /* large arguments */
    { NULL, { "0", "100" }, "osc asymp" },
    { NULL, { "1/3", "1000" }, "osc asymp" },
    { NULL, { "0", "10^4" }, "osc asymp" },
    { NULL, { "3", "10^6" }, "osc asymp" },
    { NULL, { "1/3", "10^10" }, "osc asymp" },
    { NULL, { "0", "10^15" }, "osc asymp" },
    { NULL, { "0", "10^30" }, "osc asymp" },
    { NULL, { "5/2", "10^100" }, "osc asymp special" },
    /* large orders, small arguments */
    { NULL, { "100", "1" }, "mono series" },
    { NULL, { "1000", "10" }, "mono series" },
    { NULL, { "10^4 + 1/3", "100" }, "mono series" },
    { NULL, { "10^6", "1000" }, "mono series" },
    { NULL, { "-1000 - 1/3", "10" }, "mono series" },
    /* the turning point z = nu and both sides of it */
    { NULL, { "100", "100" }, "transition turning" },
    { NULL, { "1000", "1000" }, "transition turning" },
    { NULL, { "1000 + 1/3", "1000" }, "transition turning" },
    { NULL, { "1000", "1000 + 1000^(1/3)" }, "transition turning" },
    { NULL, { "1000", "1000 - 2*1000^(1/3)" }, "transition turning" },
    { NULL, { "10^4", "10^4" }, "transition turning" },
    { NULL, { "10^5 + 1/2", "10^5" }, "transition turning" },
    { NULL, { "10^6", "10^6" }, "transition turning" },
    { NULL, { "10^6", "10^6 + 100" }, "transition turning" },
    { NULL, { "30", "100" }, "osc" },
    { NULL, { "100", "30" }, "mono" },
    { NULL, { "1000", "500" }, "mono transition" },
    { NULL, { "1000", "900" }, "mono transition" },
    { NULL, { "1000", "1100" }, "osc transition" },
    { NULL, { "1000", "2000" }, "osc transition" },
    { NULL, { "1000", "10^4" }, "osc" },
    { NULL, { "1000", "10^6" }, "osc asymp" },
    /* huge orders and arguments */
    { NULL, { "10^8", "10^8" }, "transition turning" },
    { NULL, { "10^10", "2*10^10" }, "osc transition" },
    { NULL, { "10^10", "10^10/2" }, "mono transition" },
    { NULL, { "10^15", "10^15 + 10^5" }, "transition turning" },
    { NULL, { "10^6 + 1/3", "10^9" }, "osc" },
    /* complex arguments */
    { NULL, { "0", "1 + i" }, "" },
    { NULL, { "1", "10*i" }, "mono" },
    { NULL, { "1/3", "3 - 4*i" }, "" },
    { NULL, { "0", "100 + 100*i" }, "asymp" },
    { NULL, { "1/3", "1 + 1000*i" }, "mono asymp" },
    { NULL, { "10", "10^4*(1 + i)" }, "asymp" },
    { NULL, { "100", "100*i" }, "mono transition" },
    { NULL, { "100", "100*exp(pi*i/3)" }, "transition" },
    { NULL, { "1000", "1000 + 10*i" }, "transition turning" },
    { NULL, { "1000", "1000*exp(i/10)" }, "transition turning" },
    { NULL, { "1000", "2000 - 1500*i" }, "transition" },
    { NULL, { "10^5", "10^5*(1 + i)" }, "transition" },
    /* complex orders */
    { NULL, { "i", "1" }, "" },
    { NULL, { "10*i", "10" }, "transition" },
    { NULL, { "100*i", "100" }, "transition" },
    { NULL, { "1000*i", "500" }, "transition" },
    { NULL, { "1000*i", "2000" }, "osc transition" },
    { NULL, { "10 + 10*i", "20 + 5*i" }, "transition" },
    { NULL, { "100 + 100*i", "50" }, "transition" },
    { NULL, { "1/2 + 1000*i", "1" }, "series" },
    { NULL, { "1/3 + 1000*i", "10^-3" }, "series" },
    { NULL, { "300 - 400*i", "300 + 400*i" }, "transition" },
    { NULL, { "10^4*i", "10^4" }, "transition" },
    /* the negative real axis */
    { NULL, { "0", "-10" }, "osc cut" },
    { NULL, { "3", "-5" }, "cut" },
    { NULL, { "1/3", "-5" }, "cut" },
    { NULL, { "1/3", "-5 + 10^-20*i" }, "cut" },
    { NULL, { "1000 + 1/3", "-1000" }, "cut transition turning" },
    { NULL, { "1/3", "-10^6" }, "cut osc asymp" },
    /* inexact input */
    { NULL, { "0", "10*(1 +/- 2^-30)" }, "osc" },
    { NULL, { "1/3*(1 +/- 1e-10)", "1000" }, "osc asymp" },
    { NULL, { "1000*(1 +/- 1e-20)", "1000*(1 +/- 1e-20)" }, "transition turning" },
    { NULL, { "2 +/- 1e-10", "3" }, "degenerate" },
    { NULL, { "0 +/- 1e-10", "1" }, "degenerate" },
    { NULL, { "100", "1 +/- 0.5" }, "wide series" },
    { NULL, { "1.5 +/- 0.5", "3 +/- 1" }, "wide" },
    { NULL, { "1000*" PR, "1000*" PR }, "transition turning" },
    { NULL, { "1000 + 1/3", "2000*(1 +/- 2^-(prec/2))" }, "osc transition" },
    { NULL, { "(3 + 4*i)*(1 +/- 1e-15)", "(10 - 5*i)*(1 +/- 1e-15)" }, "" },
    { NULL, { "1/3", "(-5 +/- 1e-10) + (0 +/- 1e-10)*i" }, "cut wide" },
    { NULL, { "10^6*(1 +/- 1e-12)", "10^6*(1 +/- 1e-12)" }, "transition turning" },
    { NULL, { "0", "10^30*(1 +/- 1e-25)" }, "osc asymp wide" },
};

/* I_nu(z) and K_nu(z) */
static const entry_struct entries_ik[] = {
    /* orders and arguments of moderate size */
    { NULL, { "0", "1" }, "" },
    { NULL, { "1", "1/2" }, "" },
    { NULL, { "0", "10" }, "mono" },
    { NULL, { "1/2", "3" }, "special" },
    { NULL, { "7/2", "10" }, "special" },
    { NULL, { "1/3", "5" }, "" },
    { NULL, { "-1/3", "5" }, "" },
    { NULL, { "-5/2", "2" }, "special" },
    { NULL, { "-3", "2" }, "" },
    { NULL, { "5", "5" }, "" },
    { NULL, { "pi", "sqrt(2)" }, "" },
    /* tiny arguments and orders */
    { NULL, { "0", "10^-30" }, "series" },
    { NULL, { "1", "10^-30" }, "series" },
    { NULL, { "1/3", "10^-100" }, "series" },
    { NULL, { "10^-20", "1" }, "degenerate" },
    { NULL, { "10^-30", "10^-30" }, "degenerate series" },
    /* orders near integers */
    { NULL, { "2 + 2^-30", "3" }, "degenerate" },
    { NULL, { "2 + 10^-15", "3" }, "degenerate" },
    { NULL, { "-3 + 2^-40", "2" }, "degenerate" },
    /* large arguments */
    { NULL, { "0", "100" }, "mono asymp" },
    { NULL, { "1/3", "1000" }, "mono asymp" },
    { NULL, { "0", "10^4" }, "mono asymp" },
    { NULL, { "3", "10^6" }, "mono asymp" },
    { NULL, { "1/3", "10^10" }, "mono asymp" },
    { NULL, { "0", "10^15" }, "mono asymp" },
    /* large orders, small arguments */
    { NULL, { "100", "1" }, "mono series" },
    { NULL, { "1000", "10" }, "mono series" },
    { NULL, { "10^4 + 1/3", "100" }, "mono series" },
    { NULL, { "10^6", "1000" }, "mono series" },
    { NULL, { "-1000 - 1/3", "10" }, "mono series" },
    /* argument comparable to the order */
    { NULL, { "100", "100" }, "mono transition" },
    { NULL, { "1000", "1000" }, "mono transition" },
    { NULL, { "1000 + 1/3", "1000" }, "mono transition" },
    { NULL, { "10^4", "10^4" }, "mono transition" },
    { NULL, { "10^6", "10^6" }, "mono transition" },
    { NULL, { "30", "100" }, "mono" },
    { NULL, { "100", "30" }, "mono" },
    { NULL, { "1000", "500" }, "mono transition" },
    { NULL, { "1000", "2000" }, "mono transition" },
    { NULL, { "1000", "10^4" }, "mono" },
    { NULL, { "1000", "10^6" }, "mono asymp" },
    /* huge orders and arguments */
    { NULL, { "10^8", "10^8" }, "mono transition" },
    { NULL, { "10^10", "2*10^10" }, "mono transition" },
    { NULL, { "10^15", "10^15" }, "mono transition" },
    /* complex arguments */
    { NULL, { "0", "1 + i" }, "" },
    { NULL, { "1", "10*i" }, "osc" },
    { NULL, { "1/3", "3 - 4*i" }, "" },
    { NULL, { "0", "100 + 100*i" }, "asymp" },
    { NULL, { "1/3", "1 + 1000*i" }, "osc asymp" },
    { NULL, { "10", "10^4*(1 + i)" }, "asymp" },
    { NULL, { "100", "100*i" }, "transition turning" },
    { NULL, { "1000", "1000*i" }, "transition turning" },
    { NULL, { "1000", "10 + 1000*i" }, "transition turning" },
    { NULL, { "1000", "2000*i" }, "osc transition" },
    { NULL, { "1000", "1500 + 2000*i" }, "transition" },
    { NULL, { "10^5", "10^5*(1 + i)" }, "transition" },
    /* complex orders */
    { NULL, { "i", "1" }, "" },
    { NULL, { "10*i", "10" }, "transition turning" },
    { NULL, { "100*i", "100" }, "transition turning" },
    { NULL, { "1000*i", "1000" }, "transition turning" },
    { NULL, { "1000*i", "500" }, "osc transition" },
    { NULL, { "1000*i", "2000" }, "mono transition" },
    { NULL, { "10 + 10*i", "20 + 5*i" }, "transition" },
    { NULL, { "100 + 100*i", "50" }, "transition" },
    { NULL, { "1/2 + 1000*i", "1" }, "osc series" },
    { NULL, { "1/3 + 1000*i", "10^-3" }, "osc series" },
    { NULL, { "300 - 400*i", "300 + 400*i" }, "transition" },
    { NULL, { "10^4*i", "10^4" }, "transition turning" },
    { NULL, { "10^5*i", "10^5/2" }, "osc transition" },
    /* the negative real axis */
    { NULL, { "0", "-10" }, "cut" },
    { NULL, { "3", "-5" }, "cut" },
    { NULL, { "1/3", "-5" }, "cut" },
    { NULL, { "1/3", "-5 + 10^-20*i" }, "cut" },
    { NULL, { "1000 + 1/3", "-1000" }, "cut transition" },
    /* inexact input */
    { NULL, { "0", "10*(1 +/- 2^-30)" }, "mono" },
    { NULL, { "1/3*(1 +/- 1e-10)", "1000" }, "mono asymp" },
    { NULL, { "1000*(1 +/- 1e-20)", "1000*(1 +/- 1e-20)" }, "mono transition" },
    { NULL, { "2 +/- 1e-10", "3" }, "degenerate" },
    { NULL, { "0 +/- 1e-10", "1" }, "degenerate" },
    { NULL, { "100", "1 +/- 0.5" }, "wide series" },
    { NULL, { "1.5 +/- 0.5", "3 +/- 1" }, "wide" },
    { NULL, { "1000*" PR, "1000*" PR }, "mono transition" },
    { NULL, { "1000*i*" PR, "1000*" PR }, "transition turning" },
    { NULL, { "(3 + 4*i)*(1 +/- 1e-15)", "(10 - 5*i)*(1 +/- 1e-15)" }, "" },
    { NULL, { "1/3", "(-5 +/- 1e-10) + (0 +/- 1e-10)*i" }, "cut wide" },
    { NULL, { "10^6*(1 +/- 1e-12)", "10^6*(1 +/- 1e-12)" }, "mono transition" },
};

/* Gamma(s,z) and gamma(s,z) (G, L), and the regularized Q and P where given */
static const entry_struct entries_gamma[] = {
    /* parameters and arguments of moderate size */
    { "GLQP", { "1", "1" }, "special" },
    { "GLQP", { "1/2", "1" }, "special" },
    { "G", { "0", "1" }, "degenerate" },
    { "GLQP", { "5", "3" }, "special" },
    { "GLQP", { "1/3", "2" }, "" },
    { NULL, { "-1/2", "2" }, "" },
    { "G", { "-1", "2" }, "degenerate" },
    { "G", { "-5", "3" }, "degenerate" },
    { NULL, { "pi", "sqrt(2)" }, "" },
    { "GLQP", { "100", "3" }, "series" },
    /* tiny arguments and parameters */
    { NULL, { "1/3", "10^-30" }, "series" },
    { "G", { "0", "10^-30" }, "degenerate series" },
    { "G", { "-3", "10^-10" }, "degenerate series" },
    { NULL, { "5", "10^-20" }, "series" },
    { NULL, { "10^-20", "1" }, "degenerate" },
    { NULL, { "-2 + 2^-30", "1" }, "degenerate cancel" },
    { NULL, { "2 + 2^-40", "3" }, "" },
    /* large arguments */
    { NULL, { "1/3", "100" }, "asymp" },
    { "G", { "0", "1000" }, "asymp degenerate" },
    { NULL, { "10", "10^4" }, "asymp" },
    { NULL, { "1/2", "10^6" }, "asymp special" },
    { NULL, { "1/3", "10^10" }, "asymp" },
    { NULL, { "-10 + 1/3", "10^15" }, "asymp" },
    /* large parameters, small arguments */
    { "GLQP", { "100", "1" }, "series" },
    { NULL, { "1000", "10" }, "series" },
    { NULL, { "10^4 + 1/3", "100" }, "series" },
    { NULL, { "10^6", "1" }, "series" },
    /* the transition s = z and both sides of it */
    { "GLQP", { "100", "100" }, "transition" },
    { "GLQP", { "1000", "1000" }, "transition" },
    { NULL, { "1000 + 1/3", "1000" }, "transition" },
    { NULL, { "1000", "1000 + sqrt(1000)" }, "transition" },
    { NULL, { "1000", "1000 - 3*sqrt(1000)" }, "transition" },
    { "GLQP", { "10^4", "10^4" }, "transition" },
    { "GLQP", { "10^6", "10^6" }, "transition" },
    { NULL, { "10^5", "10^5 - 1000" }, "transition" },
    { "GLQP", { "1000", "900" }, "transition" },
    { "GLQP", { "1000", "1100" }, "transition" },
    { "GLQP", { "1000", "500" }, "mono transition" },
    { "GLQP", { "1000", "2000" }, "mono transition" },
    { NULL, { "1000", "10^4" }, "mono" },
    { NULL, { "1000", "10^6" }, "asymp" },
    /* huge parameters and arguments */
    { "GLQP", { "10^8", "10^8" }, "transition" },
    { "GLQP", { "10^10", "10^10 + 10^5" }, "transition" },
    { NULL, { "10^12 + 1/3", "10^12" }, "transition" },
    { NULL, { "10^10", "10^9" }, "mono" },
    { NULL, { "10^9", "10^10" }, "mono" },
    /* negative parameters */
    { "G", { "-1000", "10" }, "degenerate" },
    { NULL, { "-1000 + 1/3", "1000" }, "transition" },
    { NULL, { "-201/2", "100" }, "transition" },
    { NULL, { "-10^4 + 1/2", "10^4" }, "transition" },
    { "G", { "-1000", "1000" }, "degenerate transition" },
    { NULL, { "-10^6 - 1/3", "10^6" }, "transition" },
    /* negative arguments */
    { NULL, { "3", "-5" }, "" },
    { NULL, { "1/3", "-5" }, "cut" },
    { NULL, { "10", "-100" }, "" },
    { NULL, { "100", "-100" }, "transition" },
    { NULL, { "1000", "-1000" }, "transition" },
    { NULL, { "1/2", "-10^4" }, "cut asymp" },
    { NULL, { "1000 + 1/3", "-1000 + 10^-20*i" }, "cut transition" },
    /* complex arguments */
    { NULL, { "1/3", "10*i" }, "osc" },
    { "G", { "0", "100*i" }, "osc asymp degenerate" },
    { NULL, { "1/2", "1000*i" }, "osc asymp special" },
    { NULL, { "1/2", "10^6*i" }, "osc asymp special" },
    { NULL, { "5", "3 + 4*i" }, "" },
    { NULL, { "100", "100*i" }, "transition" },
    { NULL, { "1000", "1000*i" }, "transition" },
    { NULL, { "1000", "1000*exp(2*pi*i/3)" }, "transition" },
    { NULL, { "1000", "-1000 + i" }, "transition" },
    { NULL, { "10", "-100 + 10^-10*i" }, "" },
    { NULL, { "1000", "500 + 800*i" }, "transition" },
    { NULL, { "10^5", "10^5*(1 + i)/sqrt(2)" }, "transition" },
    /* complex parameters */
    { NULL, { "i", "1" }, "" },
    { NULL, { "10*i", "10" }, "transition" },
    { NULL, { "1/2 + 100*i", "100" }, "osc transition" },
    { NULL, { "1/2 + 1000*i", "1000" }, "osc transition" },
    { NULL, { "1/4 + 10^4*i", "3000" }, "osc transition" },
    { NULL, { "1/4 + 10^4*i", "3000*pi*i" }, "osc transition" },
    { NULL, { "100 + 100*i", "100 - 50*i" }, "transition" },
    { NULL, { "1000*i", "1000*i" }, "transition" },
    { NULL, { "-50 + 30*i", "20 + i" }, "transition" },
    { NULL, { "1/2 + 10^5*i", "10^5" }, "osc transition" },
    /* inexact input */
    { NULL, { "1/3*(1 +/- 1e-10)", "2" }, "" },
    { "G", { "0 +/- 1e-10", "1" }, "degenerate" },
    { "G", { "-2 +/- 1e-10", "3" }, "degenerate" },
    { NULL, { "5", "3 +/- 1" }, "wide" },
    { NULL, { "1.5 +/- 0.5", "3 +/- 1" }, "wide" },
    { "GLQP", { "1000*(1 +/- 1e-20)", "1000*(1 +/- 1e-20)" }, "transition" },
    { NULL, { "10", "100*" PR }, "asymp" },
    { "GLQP", { "1000*" PR, "1000*" PR }, "transition" },
    { NULL, { "3", "(1 + i)*(1 +/- 1e-15)" }, "" },
    { NULL, { "(1/2 + 100*i)*(1 +/- 1e-15)", "100*(1 +/- 1e-15)" }, "osc transition" },
    { NULL, { "1/3", "(-5 +/- 1e-10) + (0 +/- 1e-10)*i" }, "cut wide" },
    { NULL, { "10^6*(1 +/- 1e-12)", "10^6*(1 +/- 1e-12)" }, "transition" },
};

/* 1F1(a,b,z) and U(a,b,z) */
static const entry_struct entries_kummer[] = {
    /* parameters and arguments of moderate size */
    { NULL, { "1/3", "2/3", "5" }, "" },
    { NULL, { "1", "1/2", "10" }, "" },
    { "M", { "2", "5", "-10" }, "cancel" },
    { NULL, { "1/2", "1/2", "2" }, "special" },
    { NULL, { "1/2", "3/2", "2" }, "special" },
    { NULL, { "1", "1", "1" }, "special degenerate" },
    { NULL, { "pi", "sqrt(2)", "3" }, "" },
    { NULL, { "5", "2", "3" }, "special degenerate" },
    { "U", { "2", "5", "3" }, "special poly degenerate" },
    /* polynomials */
    { NULL, { "-3", "1/2", "2" }, "poly" },
    { NULL, { "-5", "1/3", "10" }, "poly" },
    { NULL, { "-50", "1", "30" }, "poly osc cancel" },
    { NULL, { "-1000", "1", "10" }, "poly osc" },
    { NULL, { "-1000", "1", "1000" }, "poly osc transition" },
    { NULL, { "-1000", "1", "4100" }, "poly mono transition" },
    { NULL, { "-10^4", "3/2", "2*10^4" }, "poly osc transition" },
    { NULL, { "-10^5", "1/2", "10^5" }, "poly osc transition" },
    /* integer b, and b near integers */
    { "U", { "1/3", "1", "1" }, "degenerate" },
    { "U", { "1/3", "2", "1/10" }, "degenerate" },
    { "U", { "1/2", "5", "2" }, "degenerate" },
    { "U", { "31/3", "20", "5" }, "degenerate" },
    { "U", { "1/3", "-3", "2" }, "degenerate" },
    { NULL, { "1/3", "1 + 2^-30", "1" }, "degenerate cancel" },
    { NULL, { "1/3", "2 + 10^-15", "1/2" }, "degenerate cancel" },
    { NULL, { "1/3", "1 + 10^-100", "1" }, "degenerate cancel" },
    { "M", { "1", "-3 + 2^-30", "2" }, "pole" },
    { NULL, { "1", "10^-20", "2" }, "pole" },
    /* tiny arguments and parameters */
    { NULL, { "1/3", "2/3", "10^-30" }, "series" },
    { NULL, { "1", "3", "10^-10" }, "series degenerate" },
    { NULL, { "1/3", "1", "10^-30" }, "series degenerate" },
    { NULL, { "10^-20", "1", "5" }, "degenerate" },
    /* large arguments */
    { NULL, { "1/3", "2/3", "100" }, "asymp" },
    { NULL, { "1", "5/2", "10^4" }, "asymp" },
    { NULL, { "3", "7/2", "10^6" }, "asymp" },
    { NULL, { "1/3", "1", "10^10" }, "asymp degenerate" },
    { NULL, { "1/3", "2/3", "10^15" }, "asymp" },
    { "M", { "1/3", "2/3", "-100" }, "asymp cancel" },
    { "M", { "1", "5/2", "-10^4" }, "asymp cancel" },
    { "M", { "2", "7/3", "-10^6" }, "asymp cancel" },
    /* large positive parameters */
    { NULL, { "100", "200", "50" }, "mono transition" },
    { NULL, { "1000", "2000", "1000" }, "mono transition" },
    { "M", { "1000", "2000", "-1000" }, "mono transition" },
    { NULL, { "10^4", "2*10^4", "10^4" }, "mono transition" },
    { NULL, { "500", "10", "100" }, "mono transition" },
    { NULL, { "10", "500", "100" }, "mono transition" },
    { NULL, { "10", "1000 + 1/3", "10^4" }, "mono" },
    { NULL, { "1000", "1/2", "10" }, "mono" },
    { NULL, { "2000", "1000 + 1/3", "500" }, "mono transition" },
    { NULL, { "1000 + 1/3", "2000 + 1/2", "3000" }, "mono transition" },
    { NULL, { "1", "1000 + 1/2", "1000" }, "transition" },
    { NULL, { "1", "10^4 + 1/2", "10^4" }, "transition" },
    { NULL, { "500 + 5/6", "1001 + 2/3", "1000" }, "mono transition" },
    /* negative a: the oscillatory region, the turning point and beyond */
    { NULL, { "-301/3", "1/2", "50" }, "osc transition" },
    { NULL, { "-1000 - 1/2", "3", "500" }, "osc transition" },
    { NULL, { "-1000 - 1/2", "3", "3900" }, "osc transition" },
    { NULL, { "-1000 - 1/2", "3", "4008" }, "transition turning" },
    { NULL, { "-1000 - 1/2", "3", "4500" }, "mono transition" },
    { NULL, { "-10^4 - 1/3", "1", "10^4" }, "osc transition degenerate" },
    { NULL, { "-1000 - 1/3", "2000 + 1/2", "1000" }, "osc transition" },
    { NULL, { "-1000 - 1/3", "2000 + 1/2", "535" }, "transition turning" },
    { NULL, { "-1000 - 1/3", "2000 + 1/2", "300" }, "mono transition" },
    { "M", { "-21/2", "1/3", "4.6789942052908359623" }, "zero cancel" },
    /* negative b */
    { NULL, { "10", "-201/2", "20" }, "transition" },
    { NULL, { "-401/2", "-201/2", "30" }, "transition" },
    { NULL, { "100", "-401/2", "30" }, "transition" },
    /* complex arguments */
    { NULL, { "1/3", "2/3", "10*i" }, "osc" },
    { NULL, { "1", "5/2", "1000*i" }, "osc asymp" },
    { NULL, { "10", "20", "1000*(1 + i)" }, "asymp" },
    { NULL, { "50", "100", "-500 + 500*i" }, "transition" },
    { NULL, { "1000", "2000", "2000*i" }, "osc transition" },
    { NULL, { "2", "1/3", "-50 + 50*i" }, "" },
    { NULL, { "10", "1/3", "1000*exp(3*pi*i/4)" }, "asymp" },
    /* complex parameters */
    { NULL, { "i", "1", "10" }, "degenerate" },
    { NULL, { "1 + 100*i", "2", "50" }, "transition degenerate" },
    { NULL, { "1 - 10*i", "2", "40*i" }, "transition turning degenerate" },
    { NULL, { "1 - 100*i", "2", "400*i" }, "transition turning degenerate" },
    { NULL, { "1 - 100*i", "2", "100*i" }, "mono transition degenerate" },
    { NULL, { "1 - 100*i", "2", "2000*i" }, "osc transition degenerate" },
    { NULL, { "11 - 1000*i", "22", "4000*i" }, "transition degenerate" },
    { NULL, { "51 - 10*i", "102", "10*i" }, "mono transition degenerate" },
    { NULL, { "1/4 + 100*i", "1/2", "100*i" }, "transition" },
    { NULL, { "100 + 100*i", "200 - 50*i", "300 + 100*i" }, "transition" },
    { NULL, { "100 + 100*i", "50 - 20*i", "30 + 40*i" }, "transition" },
    { NULL, { "1000*i", "1 + 2000*i", "500" }, "transition" },
    { NULL, { "10*i", "1/2", "20" }, "transition" },
    { NULL, { "1/2 + 1000*i", "1", "100" }, "transition degenerate" },
    /* the branch cut of U */
    { "U", { "1", "1", "-5" }, "cut degenerate" },
    { "U", { "1/3", "1/2", "-100" }, "cut asymp" },
    { "U", { "1", "1", "-5 + 10^-10*i" }, "cut degenerate" },
    /* huge parameters and arguments */
    { NULL, { "10^6", "2*10^6", "10^6" }, "mono transition" },
    { NULL, { "1", "10^8 + 1/2", "10^8" }, "transition" },
    { NULL, { "10^10", "2*10^10 + 1/3", "10^10" }, "mono transition" },
    { NULL, { "-10^6 - 1/2", "1", "10^6" }, "osc transition degenerate" },
    { NULL, { "10^15", "10^15 + 1/2", "1" }, "mono" },
    { NULL, { "10^6", "1/2", "10^6" }, "mono transition" },
    { NULL, { "1/3", "10^8 + 1/2", "10^8" }, "transition" },
    { NULL, { "10^10 + 1/3", "1/2", "10^5" }, "mono" },
    { "M", { "1/2", "3/2", "-10^10" }, "asymp special" },
    /* inexact input */
    { NULL, { "1/3*(1 +/- 1e-10)", "2/3", "5" }, "" },
    { NULL, { "1 +/- 1e-15", "2 +/- 1e-15", "100 +/- 1e-13" }, "asymp degenerate" },
    { NULL, { "-5 +/- 1e-10", "2", "3" }, "poly degenerate" },
    { NULL, { "1/3", "1 +/- 1e-10", "1" }, "degenerate" },
    { "M", { "1", "-2 +/- 1e-10", "1" }, "pole" },
    { NULL, { "1.5 +/- 0.5", "3 +/- 0.5", "2 +/- 1" }, "wide" },
    { NULL, { "1000*" PR, "2000*" PR, "1000*" PR }, "mono transition" },
    { NULL, { "(-1000 - 1/2)*" PR, "3*" PR, "500*" PR }, "osc transition" },
    { NULL, { "(-1000 - 1/2)*(1 +/- 1e-20)", "3", "4008*(1 +/- 1e-20)" }, "transition turning" },
    { NULL, { "(1 - 100*i)*(1 +/- 1e-15)", "2", "400*i*(1 +/- 1e-15)" }, "transition turning degenerate" },
    { NULL, { "(100 + 100*i)*" PR, "(200 - 50*i)*" PR, "(300 + 100*i)*" PR }, "transition" },
    { "U", { "1", "1", "(-5 +/- 1e-10) + (0 +/- 1e-10)*i" }, "cut wide degenerate" },
    { NULL, { "10^6*(1 +/- 1e-12)", "2*10^6*(1 +/- 1e-12)", "10^6*(1 +/- 1e-12)" }, "mono transition" },
    /* polynomials with exact parameters and a rounded argument (the
       Taylor expansion in z) */
    { NULL, { "-50", "1", "100/3" }, "poly osc cancel" },
    { NULL, { "-1000", "1", "31/3" }, "poly osc cancel" },
    { NULL, { "-300", "3/8", "151/3" }, "poly osc cancel" },
    { "M", { "301", "1", "-100/3" }, "poly cancel" },
    { "U", { "1000", "2000", "1000/3" }, "poly cancel" },
};

/* 2F1(a,b,c,z) */
static const entry_struct entries_gauss[] = {
    /* parameters of moderate size, z inside the unit disc */
    { NULL, { "1", "1", "2", "1/2" }, "special" },
    { NULL, { "1/2", "1/2", "1", "1/4" }, "special degenerate" },
    { NULL, { "1/3", "2/3", "5/6", "1/3" }, "" },
    { NULL, { "1", "2", "3", "-1" }, "special" },
    { NULL, { "1/2", "1/3", "1/4", "9/10" }, "" },
    { NULL, { "2", "3", "7/2", "99/100" }, "" },
    { NULL, { "pi", "sqrt(2)", "exp(1)", "1/3" }, "" },
    { NULL, { "1/3", "1/4", "1/5", "-3" }, "" },
    { NULL, { "1/3", "2/3", "5/6", "10^-30" }, "series" },
    /* large arguments */
    { NULL, { "1", "1", "2", "-10" }, "special degenerate" },
    { NULL, { "1/3", "1/4", "1/5", "-1000" }, "asymp" },
    { NULL, { "1", "1", "2", "10^6" }, "cut asymp degenerate" },
    { NULL, { "1/3", "2/3", "1", "-10^10" }, "asymp" },
    { NULL, { "1/3", "1/4", "1/5", "10^10*i" }, "asymp" },
    /* z near 1 */
    { NULL, { "1/3", "2/3", "1", "1 - 10^-20" }, "degenerate" },
    { NULL, { "1/3", "1/4", "2", "1" }, "" },
    { NULL, { "1/3", "1/4", "1/5", "1 - 10^-10" }, "" },
    { NULL, { "1/3", "1/4", "2", "1 + 10^-10" }, "cut" },
    /* z near exp(+-pi i/3) */
    { NULL, { "1/3", "2/3", "1", "exp(pi*i/3)" }, "corner degenerate" },
    { NULL, { "1/3", "1/4", "1/5", "(1 + sqrt(3)*i)/2" }, "corner" },
    { NULL, { "1/3", "1/4", "1/5", "(1 + sqrt(3)*i)/2*(1 + 10^-3)" }, "corner" },
    { NULL, { "1/3", "1/4", "1/5", "1/2 - 7/8*i" }, "corner" },
    /* the unit circle and the branch cut */
    { NULL, { "1/3", "1/4", "2", "i" }, "" },
    { NULL, { "1/3", "1/4", "2", "exp(3*i)" }, "" },
    { NULL, { "1/3", "1/4", "2", "-1" }, "" },
    { NULL, { "1/3", "1/4", "1/5", "2" }, "cut" },
    { NULL, { "1/3", "1/4", "1/5", "2 + 10^-10*i" }, "cut" },
    { NULL, { "1/3", "1/4", "1/5", "3/2 - 10^-30*i" }, "cut" },
    /* integer differences of the parameters, and near them */
    { NULL, { "1/3", "4/3", "1/5", "-5" }, "degenerate" },
    { NULL, { "2", "5", "1/3", "10*i" }, "degenerate" },
    { NULL, { "1/3", "2/3", "2", "9/10" }, "degenerate" },
    { NULL, { "1/3", "1/3 + 10^-15", "1/5", "-5" }, "degenerate cancel" },
    { NULL, { "1/3", "2/3", "1 + 10^-20", "19/20" }, "degenerate cancel" },
    { NULL, { "1/3", "2/3", "1 + 2^-40", "1/2 + i" }, "degenerate cancel corner" },
    { NULL, { "1", "1", "-3 + 2^-30", "1/3" }, "pole" },
    /* polynomials */
    { NULL, { "-10", "1/3", "1/5", "2" }, "poly" },
    { NULL, { "-20", "1/3", "1/5", "-1000" }, "poly" },
    { NULL, { "-50", "51", "1", "1/3" }, "poly osc" },
    { NULL, { "-100", "101", "1", "3" }, "poly mono cancel" },
    { NULL, { "-100", "101", "1", "-2" }, "poly mono" },
    { NULL, { "-100", "1/2", "-200", "1/3" }, "poly" },
    { NULL, { "-200", "200 + 7/3", "3/2", "3/10" }, "poly osc" },
    { NULL, { "-1000", "1001", "1", "1/4" }, "poly osc" },
    { NULL, { "-1000", "1001", "1", "10^-4" }, "poly osc" },
    { NULL, { "-10^4", "10^4 + 1", "1", "1/2" }, "poly osc" },
    { NULL, { "-10^5", "10^5 + 1", "1", "1/3" }, "poly osc" },
    /* Legendre and Jacobi functions of large degree (a < 0 < c < b) */
    { NULL, { "-201/2", "203/2", "1", "1/4" }, "osc" },
    { NULL, { "-1000 - 1/2", "1001 + 1/2", "1", "1/4" }, "osc" },
    { NULL, { "-1000 - 1/2", "1001 + 1/2", "1", "-1/4" }, "mono" },
    { NULL, { "-10^4 - 1/2", "10^4 + 3/2", "1", "1/3" }, "osc" },
    { NULL, { "1000", "-1000 + 1/2", "1/2", "3/10" }, "osc" },
    { NULL, { "-1000 - 1/3", "1000 + 1/2", "1/3", "9/10" }, "osc" },
    { NULL, { "-21/2", "23/2", "1", "0.67455270640740607053" }, "zero cancel" },
    /* large parameters with c > b > 0 */
    { NULL, { "10", "20", "50", "1/2" }, "mono" },
    { NULL, { "100", "200", "500", "1/2" }, "mono" },
    { NULL, { "1000", "2000", "5000", "1/2" }, "mono" },
    { NULL, { "1000", "2000", "5000", "9/10" }, "mono" },
    { NULL, { "1000", "2000", "5000", "-3" }, "mono" },
    { NULL, { "1000", "2000", "5000", "-20" }, "mono" },
    { NULL, { "1000", "2000", "5000", "6/5" }, "cut turning" },
    { NULL, { "1000", "2000", "5000", "5/2" }, "cut osc" },
    { NULL, { "1000", "2000", "5000", "21" }, "cut turning" },
    { NULL, { "1000", "2000", "5000", "5 + i/10" }, "osc" },
    { NULL, { "1000", "2000", "5000", "2 + 3*i" }, "" },
    { NULL, { "1000", "10", "3000", "9/10" }, "mono" },
    { NULL, { "1", "1000", "2000", "99/100" }, "mono" },
    { NULL, { "5", "2", "400", "11/10*i" }, "" },
    { NULL, { "1000", "1/3", "2000", "3/4 + i" }, "" },
    /* large c, large a and b, c near a + b */
    { NULL, { "1/3", "2/3", "1000", "1/2" }, "series" },
    { NULL, { "1", "1", "10^6", "-10^6" }, "transition" },
    { NULL, { "1/3", "2/3", "10^6 + 1/2", "10^5*i" }, "transition" },
    { NULL, { "100", "100", "1", "1/10" }, "mono degenerate" },
    { NULL, { "1000", "1001", "2", "1/100" }, "mono degenerate" },
    { NULL, { "100 + 1/3", "100 + 1/4", "1/5", "-1/2" }, "mono" },
    { NULL, { "1000", "1000 + 1/3", "1/2", "-10" }, "mono" },
    { NULL, { "1000", "1000", "2000 + 1/3", "99/100" }, "mono" },
    /* complex parameters */
    { NULL, { "1/2 + 100*i", "1/2 - 100*i", "1", "-3" }, "osc degenerate" },
    { NULL, { "1/2 + 1000*i", "1/2 - 1000*i", "1", "-1/2" }, "osc degenerate" },
    { NULL, { "1/2 + 100*i", "1/2 - 100*i", "1", "3/4" }, "mono degenerate" },
    { NULL, { "100 + 100*i", "50", "300 - 30*i", "1/3 + i/2" }, "" },
    { NULL, { "10*i", "20*i", "30*i", "3/10" }, "" },
    { NULL, { "10*i", "1/2", "1/2 + 20*i", "-1 + i" }, "" },
    { NULL, { "1000*i", "1/3", "1/2 + 2000*i", "1/2" }, "" },
    { NULL, { "1/3", "1/4", "1000*i", "500*i" }, "transition" },
    /* huge parameters */
    { NULL, { "10^5", "2*10^5", "5*10^5", "1/2" }, "mono" },
    { NULL, { "10^6", "1/3", "2*10^6 + 1/2", "-1" }, "mono" },
    { NULL, { "1", "10^8", "10^8 + 10^4", "9/10" }, "mono" },
    { NULL, { "1/3", "2/3", "10^10", "10^9*i" }, "transition" },
    { NULL, { "-10^6", "10^6 + 1", "1", "1/2" }, "poly osc" },
    /* inexact input */
    { NULL, { "1/3*(1 +/- 1e-10)", "2/3", "5/6", "1/3" }, "" },
    { NULL, { "1 +/- 1e-15", "2 +/- 1e-15", "3 +/- 1e-15", "(1/2)*(1 +/- 1e-15)" }, "degenerate" },
    { NULL, { "1/3", "4/3 +/- 1e-10", "1/5", "-5" }, "degenerate" },
    { NULL, { "1/3", "2/3", "1 +/- 1e-12", "19/20" }, "degenerate" },
    { NULL, { "-10 +/- 1e-10", "1/3", "1/5", "1/2" }, "poly degenerate" },
    { NULL, { "1", "1", "-3 +/- 1e-10", "1/3" }, "pole" },
    { NULL, { "1.5 +/- 0.5", "2 +/- 0.5", "5 +/- 0.5", "0.3 +/- 0.2" }, "wide" },
    { NULL, { "1000*" PR, "2000*" PR, "5000*" PR, "(1/2)*" PR }, "mono" },
    { NULL, { "(-1000 - 1/2)*" PR, "(1001 + 1/2)*" PR, "1", "(1/4)*" PR }, "osc" },
    { NULL, { "1000*(1 +/- 1e-20)", "2000*(1 +/- 1e-20)", "5000*(1 +/- 1e-20)", "-3*(1 +/- 1e-20)" }, "mono" },
    { NULL, { "(100 + 100*i)*(1 +/- 1e-15)", "50", "(300 - 30*i)*(1 +/- 1e-15)", "(1/3 + i/2)*(1 +/- 1e-15)" }, "" },
    { NULL, { "1/3", "1/4", "1/5", "(2 +/- 1e-10) + (0 +/- 1e-10)*i" }, "cut wide" },
    { NULL, { "1/3", "2/3", "1", "1 +/- 1e-10" }, "pole wide" },
    { NULL, { "10^5*(1 +/- 1e-12)", "2*10^5*(1 +/- 1e-12)", "5*10^5*(1 +/- 1e-12)", "1/2" }, "mono" },
    /* a, b > c > 0: polynomials in z/(z - 1) if c - a or c - b is an
       integer, which cancel for z < 0, where the function oscillates */
    { NULL, { "423", "412", "153", "1/2" }, "poly mono" },
    { NULL, { "423", "412", "153", "-1" }, "poly osc cancel" },
    { NULL, { "423", "412", "153", "-1/10" }, "poly osc cancel" },
    { NULL, { "423", "412", "153", "i" }, "poly cancel" },
    { NULL, { "423 + 1/3", "412 + 1/7", "153 + 1/5", "-1" }, "osc cancel" },
    { NULL, { "423 + 10*i", "412", "153 - 5*i", "-1 + i/2" }, "osc cancel" },
    { NULL, { "42300 + 1/3", "41200", "15300", "-3/2" }, "osc cancel" },
    { NULL, { "423*(1 +/- 2^-53)", "412*(1 +/- 2^-53)", "153*(1 +/- 2^-53)", "-1*(1 +/- 2^-53)" }, "osc cancel" },
    { NULL, { "423*" PR, "412*" PR, "153*" PR, "-1*" PR }, "osc cancel" },
    /* polynomials with exact parameters and a rounded argument (the
       Taylor expansion in z) */
    { NULL, { "-1000", "1001", "1", "1/3" }, "poly osc cancel" },
    { NULL, { "-300", "3/8", "5/8", "1/5" }, "poly cancel" },
    { NULL, { "423", "412", "153", "11/10" }, "poly cancel" },
    { NULL, { "37", "12", "48", "(1694109845964387/2^51) + (-741719999514547/2^50)*i" }, "corner" },
    /* the choice among the series, the connection formulas and the
       corner algorithm by their estimated costs: the other examples of
       issues 2964 and 2967; limit cases (c - a - b an integer) at
       |1 - z| = 0.45 and 0.88, the second one with complex parameters;
       z near 1 with parameters for which the series at z has few
       terms; a connection formula with twice the terms of the two
       others */
    { NULL, { "5", "7", "33", "(99/100)*i" }, "transition" },
    { NULL, { "32", "14", "48", "(6573620216893439/2^53) + (3078845386612347/2^52)*i" }, "corner" },
    { NULL, { "423", "412", "153", "(-1635325596928791/2^46) + (7558374922216237/2^49)*i" }, "cancel" },
    { NULL, { "5", "5", "12", "(729/1024) + (349/1024)*i" }, "transition" },
    { NULL, { "1/2 + i/4", "5/4", "15/4 + i/4", "3/8 + 5/8*i" }, "transition" },
    { NULL, { "-251/32", "11/64", "3899/256", "(251/256) + (-3/256)*i" }, "transition" },
    { NULL, { "849/32", "-1069/64", "989/128", "(357/16) + (73/16)*i" }, "cancel" },
};

typedef struct
{
    const char * name;
    const char * ids;           /* the functions of the table */
    const char * default_sel;
    const entry_struct * entries;
    slong num;
}
group_struct;

#define TABLE(t) t, (slong) (sizeof(t) / sizeof(t[0]))

static const group_struct groups[] = {
    { "J, Y", "JY", "JY", TABLE(entries_jy) },
    { "I, K", "IK", "IK", TABLE(entries_ik) },
    { "incomplete gamma", "GLQP", "GL", TABLE(entries_gamma) },
    { "1F1, U", "MU", "MU", TABLE(entries_kummer) },
    { "2F1", "F", "F", TABLE(entries_gauss) },
};

#define NUM_GROUPS ((slong) (sizeof(groups) / sizeof(groups[0])))

/* ------------------------------------------------------------------------- */
/*  Options                                                                  */
/* ------------------------------------------------------------------------- */

typedef struct
{
    slong precs[MAX_PRECS];
    int nprecs;
    const char * func_sel;
    const char * class_sel;
    const char * case_sel;
    const char * inputs[MAX_INPUTS];
    int ninputs;
    const char * alg_sel;
    int list;
    int list_algs;
    int ref;
    int values;
    const char * csv_name;
    FILE * csv;
    int summary_only;
    int use_fork;
    int verbose;
    int threads;
    double mintime;
    int warm;
    double timeout;
    double ratio;
    slong slack;
    slong maxmem;
    slong maxwp;
}
options_struct;

static options_struct opt;

static void
usage(void)
{
    printf(
"Usage: p-catalog [options]\n"
"\n"
"Times the functions of acb_hypgeom on a catalog of inputs and prints the\n"
"accuracy of the results.\n"
"\n"
"Selection:\n"
"  -func LIST     functions: names or letters (see -algs), or the groups\n"
"                 bessel, gamma, hypgeom, all (default all)\n"
"  -prec LIST     precisions in bits (default 64,333,3333)\n"
"  -class EXPR    classes of inputs (see -list): a+b for both, a,b for\n"
"                 either, ^a for not a; for example  turning+^inexact,huge\n"
"  -case LIST     cases by identifier (J12,K3) or all cases of a letter (F)\n"
"  -input CALL    an input which is not in the catalog, in place of it, for\n"
"                 example  -input \"2f1(1000, 1/3, 2000 +/- 1e-10, 3/4 + i)\";\n"
"                 the option can be repeated (cases X1, X2, ...)\n"
"  -alg LIST      algorithms to compare with the default function: names\n"
"                 (see -algs), all, loop:NAME, loop (= loop:default) or\n"
"                 loop:all.  A loop calls the algorithm at the precisions\n"
"                 prec, 2 prec, 4 prec, ... with the inputs evaluated at\n"
"                 each, until the result is accurate to prec - slack bits\n"
"\n"
"Output:\n"
"  -list          print the selected cases with their classes and exit\n"
"  -algs          print the functions and their algorithms and exit\n"
"  -values        print the values of the inputs and of the results\n"
"  -ref           compute reference values at the midpoints of the inputs and\n"
"                 at points on their boundaries (with the default function at\n"
"                 higher precision): checks that the results contain these\n"
"                 values, and prints the accuracy which the inputs allow\n"
"  -summary       print only the summary\n"
"  -csv FILE      also write one line of comma-separated values per\n"
"                 measurement to FILE (times in seconds)\n"
"\n"
"Measurement:\n"
"  -mintime S     repeat each call for at least S seconds (default 0.01)\n"
"  -warm          also repeat a call which takes more than half a second,\n"
"                 once: its time with the caches filled, like the others\n"
"                 (without this, the time of the first call in a new\n"
"                 process, with the quadrature nodes of an integral)\n"
"  -timeout S     time limit for one measurement (default 10); an algorithm\n"
"                 other than the default is also stopped after 100 times the\n"
"                 time of the default function (at least 0.1 seconds)\n"
"  -maxmem MB     memory limit for one measurement (default 4096; 0 for none)\n"
"  -slack B       a result is accurate if it has prec - B bits (default 10),\n"
"                 or B bits less than the inputs allow with -ref\n"
"  -maxwp N       loops stop at the precision N prec (default 32)\n"
"  -ratio R       list an algorithm or loop as faster than the default\n"
"                 function from this ratio of times (default 1.5)\n"
"  -threads N     number of threads for FLINT (default 1)\n"
"  -nofork        measure in this process: no time or memory limit (for use\n"
"                 with a profiler or a debugger)\n"
"  -verbose       do not discard the error output of the measurements\n"
"\n"
"Each line gives the time of one call (of the whole loop for a loop), with\n"
"the caches of FLINT filled (the first call in a new process is noted where\n"
"it takes more than five times as long and more than a millisecond), and the\n"
"relative accuracy of the result in bits; with -ref, the line of the default\n"
"function also gives the accuracy which the inputs allow, and with -alg the\n"
"other lines give the time over the time of the default function.  All\n"
"results of a case, at all precisions and from all algorithms, are checked\n"
"to overlap.  The summary counts the accurate results of the default\n"
"functions and lists the cases in which they are not accurate and, with\n"
"-alg, those in which another algorithm or a loop gives an accurate result\n"
"faster.\n");
}

/* ------------------------------------------------------------------------- */
/*  Utilities                                                                */
/* ------------------------------------------------------------------------- */

static double
now(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + 1e-6 * tv.tv_usec;
}

static int
is_ident_char(char c)
{
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_';
}

/* Evaluates an expression at the given precision, with "prec" replaced
   by the precision.  Returns 0 on failure. */
static int
eval_expr(acb_t res, const char * expr, slong prec)
{
    gr_ctx_t ctx;
    char * s;
    char num[32];
    size_t len, numlen, i, j;
    int status;

    len = strlen(expr);
    snprintf(num, sizeof(num), "%ld", (long) prec);
    numlen = strlen(num);
    s = flint_malloc(len * (numlen + 1) + 1);

    for (i = 0, j = 0; i < len; )
    {
        if (strncmp(expr + i, "prec", 4) == 0
            && (i == 0 || !is_ident_char(expr[i - 1])) && !is_ident_char(expr[i + 4]))
        {
            memcpy(s + j, num, numlen);
            j += numlen;
            i += 4;
        }
        else
        {
            s[j++] = expr[i++];
        }
    }
    s[j] = '\0';

    gr_ctx_init_complex_acb(ctx, prec);
    status = gr_set_str(res, s, ctx);
    gr_ctx_clear(ctx);
    flint_free(s);

    return status == GR_SUCCESS && acb_is_finite(res);
}

/* A growing buffer of bytes for the results of a measurement */

typedef struct
{
    char * data;
    size_t len;
    size_t alloc;
    size_t pos;
}
buf_struct;

typedef buf_struct buf_t[1];

static void
buf_init(buf_t b)
{
    b->data = NULL;
    b->len = b->alloc = b->pos = 0;
}

static void
buf_clear(buf_t b)
{
    flint_free(b->data);
}

static void
buf_put(buf_t b, const void * src, size_t n)
{
    if (b->len + n > b->alloc)
    {
        b->alloc = 2 * (b->len + n) + 64;
        b->data = flint_realloc(b->data, b->alloc);
    }
    memcpy(b->data + b->len, src, n);
    b->len += n;
}

static int
buf_get(buf_t b, void * dest, size_t n)
{
    if (b->pos + n > b->len)
        return 0;
    memcpy(dest, b->data + b->pos, n);
    b->pos += n;
    return 1;
}

static void
buf_put_arb(buf_t b, const arb_t x)
{
    char * s = arb_dump_str(x);
    size_t n = strlen(s) + 1;
    buf_put(b, &n, sizeof(n));
    buf_put(b, s, n);
    flint_free(s);
}

static int
buf_get_arb(buf_t b, arb_t x)
{
    size_t n;
    if (!buf_get(b, &n, sizeof(n)) || b->pos + n > b->len || n == 0
        || b->data[b->pos + n - 1] != '\0')
        return 0;
    if (arb_load_str(x, b->data + b->pos))
        return 0;
    b->pos += n;
    return 1;
}

static void
buf_put_acb(buf_t b, const acb_t x)
{
    buf_put_arb(b, acb_realref(x));
    buf_put_arb(b, acb_imagref(x));
}

static int
buf_get_acb(buf_t b, acb_t x)
{
    return buf_get_arb(b, acb_realref(x)) && buf_get_arb(b, acb_imagref(x));
}

/* Runs work(out, arg), in a child process with a time limit if possible.
   Returns 0 on success, 1 on timeout and 2 if the child did not finish
   (with the signal in *sig, or 0). */

#define RUN_OK 0
#define RUN_TIMEOUT 1
#define RUN_CRASH 2

typedef void (*work_func)(buf_t out, void * arg);

static int
run_isolated(buf_t out, work_func work, void * arg, double timeout, int * sig)
{
    *sig = 0;

#if CATALOG_HAVE_FORK
    if (opt.use_fork)
    {
        int fd[2], status, result;
        pid_t pid;
        double deadline;
        char tmp[4096];

        fflush(stdout);
        fflush(stderr);

        if (pipe(fd) != 0)
        {
            perror("pipe");
            exit(1);
        }

        pid = fork();

        if (pid < 0)
        {
            perror("fork");
            exit(1);
        }

        if (pid == 0)
        {
            buf_t b;
            size_t done;

            close(fd[0]);

            /* (in case the parent is killed) */
            alarm((unsigned int) FLINT_MIN(timeout + 5.0, 1e6));

            if (opt.maxmem > 0)
            {
                struct rlimit lim;
                lim.rlim_cur = lim.rlim_max = (rlim_t) opt.maxmem * 1024 * 1024;
                setrlimit(RLIMIT_AS, &lim);
            }

            if (!opt.verbose)
            {
                int nul = open("/dev/null", O_WRONLY);
                if (nul >= 0)
                {
                    dup2(nul, 2);
                    close(nul);
                }
            }

            buf_init(b);
            work(b, arg);

            for (done = 0; done < b->len; )
            {
                ssize_t n = write(fd[1], b->data + done, b->len - done);
                if (n < 0)
                {
                    if (errno == EINTR)
                        continue;
                    _exit(3);
                }
                done += n;
            }

            close(fd[1]);
            _exit(0);
        }

        close(fd[1]);
        deadline = now() + timeout;
        result = RUN_OK;

        while (1)
        {
            struct pollfd pfd;
            double left = deadline - now();
            int r;

            if (left <= 0)
            {
                result = RUN_TIMEOUT;
                break;
            }

            pfd.fd = fd[0];
            pfd.events = POLLIN;
            r = poll(&pfd, 1, (int) FLINT_MIN(left * 1000 + 1, 1e6));

            if (r < 0 && errno != EINTR)
            {
                result = RUN_CRASH;
                break;
            }

            if (r > 0)
            {
                ssize_t n = read(fd[0], tmp, sizeof(tmp));
                if (n > 0)
                    buf_put(out, tmp, n);
                else if (n == 0)
                    break;
                else if (errno != EINTR)
                {
                    result = RUN_CRASH;
                    break;
                }
            }
        }

        close(fd[0]);

        if (result != RUN_OK)
            kill(pid, SIGKILL);

        while (waitpid(pid, &status, 0) < 0 && errno == EINTR)
            ;

        if (result == RUN_OK)
        {
            if (WIFSIGNALED(status))
            {
                *sig = WTERMSIG(status);
                result = RUN_CRASH;
            }
            else if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
            {
                result = RUN_CRASH;
            }
        }

        return result;
    }
#endif

    (void) timeout;
    work(out, arg);
    return RUN_OK;
}

/* ------------------------------------------------------------------------- */
/*  Cases                                                                    */
/* ------------------------------------------------------------------------- */

typedef struct
{
    const func_struct * func;
    const entry_struct * entry;
    char id[32];
    char tags[256];
}
case_struct;

static case_struct * cases;
static slong num_cases;

static const func_struct *
func_from_id(char id)
{
    slong i;
    for (i = 0; i < NUM_FUNCS; i++)
        if (funcs[i].id == id)
            return funcs + i;
    return NULL;
}

static void
add_tag(char * tags, const char * tag)
{
    if (tags[0] != '\0')
        strcat(tags, " ");
    strcat(tags, tag);
}

static int
has_tag(const char * tags, const char * tag, size_t len)
{
    const char * s = tags;

    while (*s != '\0')
    {
        const char * e = s;
        while (*e != '\0' && *e != ' ')
            e++;
        if ((size_t) (e - s) == len && strncmp(s, tag, len) == 0)
            return 1;
        s = e;
        while (*s == ' ')
            s++;
    }

    return 0;
}

static double
acb_log2_abs(const acb_t x)
{
    mag_t m;
    double r;
    mag_init(m);
    acb_get_mag(m, x);
    r = mag_is_zero(m) ? -1e300 : mag_get_d_log2_approx(m);
    mag_clear(m);
    return r;
}

/* The classes which are found from the inputs */
static int
classify(char * tags, const func_struct * func, const entry_struct * entry)
{
    acb_struct x[MAX_ARGS];
    acb_t t;
    slong i;
    int n = func->nargs, ok = 1;
    int cplx = 0, exact = 1, radius = 0, precrad = 0;
    int isint = 0, halfint = 0, nearint = 0, ptiny = 0;
    double pmax = -1e300, zmag;

    tags[0] = '\0';
    acb_init(t);

    for (i = 0; i < n; i++)
    {
        acb_init(x + i);
        if (!eval_expr(x + i, entry->args[i], 128))
            ok = 0;
        if (strstr(entry->args[i], "+/-") != NULL)
        {
            radius = 1;
            if (strstr(entry->args[i], "prec") != NULL)
                precrad = 1;
        }
    }

    if (ok)
    {
        for (i = 0; i < n; i++)
        {
            if (!acb_is_real(x + i))
                cplx = 1;
            if (!acb_is_exact(x + i))
                exact = 0;
        }

        /* An expression which is not exact at 128 bits can be exact at a
           higher precision */
        if (!exact && !radius)
        {
            exact = 1;
            for (i = 0; i < n; i++)
                if (!eval_expr(t, entry->args[i], 4096) || !acb_is_exact(t))
                    exact = 0;
        }

        for (i = 0; i < n - 1; i++)
        {
            double mag = acb_log2_abs(x + i);

            pmax = FLINT_MAX(pmax, mag);
            if (mag > -1e299 && mag < -10)
                ptiny = 1;

            if (acb_is_int(x + i))
                isint = 1;
            else
            {
                acb_mul_2exp_si(t, x + i, 1);
                if (acb_is_int(t))
                    halfint = 1;
                else if (acb_contains_int(x + i))
                    nearint = 1;
                else
                {
                    /* distance from the midpoint to the integers */
                    acb_zero(t);
                    arf_floor(arb_midref(acb_realref(t)), arb_midref(acb_realref(x + i)));
                    acb_sub(t, x + i, t, 128);
                    if (arf_cmp_d(arb_midref(acb_realref(t)), 0.5) > 0)
                        arb_sub_ui(acb_realref(t), acb_realref(t), 1, 128);
                    if (acb_log2_abs(t) < -10)
                        nearint = 1;
                }
            }
        }

        zmag = acb_log2_abs(x + n - 1);

        add_tag(tags, radius ? "inexact" : (exact ? "exact" : "rounded"));
        if (radius)
            add_tag(tags, precrad ? "precrad" : "fixedrad");
        add_tag(tags, cplx ? "complex" : "real");
        if (isint)
            add_tag(tags, "int");
        if (halfint)
            add_tag(tags, "halfint");
        if (nearint)
            add_tag(tags, "nearint");
        if (pmax > -1e299 && pmax < -10)
            add_tag(tags, "ptiny");
        else
        {
            if (ptiny)
                add_tag(tags, "ptiny");
            add_tag(tags, pmax < 6 ? "pnormal" : (pmax < log(30000.0) / log(2.0) ? "plarge" : "phuge"));
        }
        if (zmag < -1e299)
            add_tag(tags, "zzero");
        else
            add_tag(tags, zmag < -10 ? "ztiny" : (zmag < 6 ? "znormal" :
                (zmag < log(30000.0) / log(2.0) ? "zlarge" : "zhuge")));
    }
    else
    {
        add_tag(tags, "INVALID");
    }

    if (entry->tags[0] != '\0')
        add_tag(tags, entry->tags);

    for (i = 0; i < n; i++)
        acb_clear(x + i);
    acb_clear(t);

    return ok;
}

static int
build_cases(void)
{
    slong g, e, alloc = 0;
    int ok = 1;

    cases = NULL;
    num_cases = 0;

    for (g = 0; g < NUM_GROUPS; g++)
    {
        const char * c;

        for (c = groups[g].ids; *c != '\0'; c++)
        {
            for (e = 0; e < groups[g].num; e++)
            {
                const entry_struct * entry = groups[g].entries + e;
                const char * sel = entry->sel != NULL ? entry->sel : groups[g].default_sel;
                case_struct * cs;

                if (strchr(sel, *c) == NULL)
                    continue;

                if (num_cases == alloc)
                {
                    alloc = 2 * alloc + 64;
                    cases = flint_realloc(cases, alloc * sizeof(case_struct));
                }

                cs = cases + num_cases;
                cs->func = func_from_id(*c);
                cs->entry = entry;
                snprintf(cs->id, sizeof(cs->id), "%c%d", *c, (int) (e + 1));

                if (!classify(cs->tags, cs->func, entry))
                {
                    printf("error: the inputs of case %s cannot be evaluated\n", cs->id);
                    ok = 0;
                }

                num_cases++;
            }
        }
    }

    return ok;
}

/* The cases of the option -input, "name(arg, ..., arg)", in place of the
   catalog */
static int
build_inputs(void)
{
    slong i, f;
    int k;

    cases = flint_malloc(opt.ninputs * sizeof(case_struct));
    num_cases = 0;

    for (i = 0; i < opt.ninputs; i++)
    {
        const char * str = opt.inputs[i];
        const char * open = strchr(str, '(');
        const char * end = str + strlen(str);
        const char * a, * b;
        const func_struct * func = NULL;
        entry_struct * entry;
        case_struct * cs = cases + num_cases;
        size_t len;
        int depth, nargs = 0;

        while (end > str && end[-1] == ' ')
            end--;

        len = (open == NULL) ? 0 : (size_t) (open - str);
        while (len > 0 && str[len - 1] == ' ')
            len--;

        for (f = 0; f < NUM_FUNCS && open != NULL; f++)
            if ((strlen(funcs[f].name) == len && strncmp(funcs[f].name, str, len) == 0)
                || (len == 1 && str[0] == funcs[f].id))
                func = funcs + f;

        if (func == NULL || end[-1] != ')')
        {
            printf("error: not of the form function(arguments), with a function of -algs: %s\n", str);
            return 0;
        }

        /* the arguments, separated by the commas outside of parentheses */
        entry = flint_calloc(1, sizeof(entry_struct));
        entry->tags = "";

        for (a = open + 1, depth = 0; a < end; a = b + 1)
        {
            char * arg;

            for (b = a; b < end - 1 && !(*b == ',' && depth == 0); b++)
                depth += (*b == '(') - (*b == ')');

            if (nargs == func->nargs)
            {
                nargs++;
                break;
            }

            arg = flint_malloc(b - a + 1);
            memcpy(arg, a, b - a);
            arg[b - a] = '\0';
            for (k = 0; arg[k] == ' '; k++) ;
            memmove(arg, arg + k, strlen(arg + k) + 1);
            for (k = strlen(arg); k > 0 && arg[k - 1] == ' '; k--)
                arg[k - 1] = '\0';
            entry->args[nargs++] = arg;
        }

        if (nargs != func->nargs)
        {
            printf("error: %s takes %d arguments: %s\n", func->name, func->nargs, str);
            return 0;
        }

        cs->func = func;
        cs->entry = entry;
        snprintf(cs->id, sizeof(cs->id), "X%d", (int) (i + 1));

        if (!classify(cs->tags, func, entry))
        {
            printf("error: the arguments cannot be evaluated: %s\n", str);
            return 0;
        }

        num_cases++;
    }

    return 1;
}

/* Whether name is one of the comma-separated items of list */
static int
in_list(const char * list, const char * name)
{
    size_t len = strlen(name);
    const char * s = list;

    while (*s != '\0')
    {
        const char * e = strchr(s, ',');
        size_t n = (e == NULL) ? strlen(s) : (size_t) (e - s);
        if (n == len && strncmp(s, name, len) == 0)
            return 1;
        if (e == NULL)
            break;
        s = e + 1;
    }

    return 0;
}

static int
func_selected(const func_struct * func)
{
    char id[2];
    const char * sel = opt.func_sel;

    if (sel == NULL || in_list(sel, "all"))
        return 1;

    id[0] = func->id;
    id[1] = '\0';

    if (in_list(sel, func->name) || in_list(sel, id))
        return 1;
    if (in_list(sel, "bessel") && strchr("JYIK", func->id) != NULL)
        return 1;
    if (in_list(sel, "gamma") && strchr("GLQP", func->id) != NULL)
        return 1;
    if (in_list(sel, "hypgeom") && strchr("MUF", func->id) != NULL)
        return 1;
    if (func->id == 'M' && (in_list(sel, "m") || in_list(sel, "1F1")))
        return 1;
    if (func->id == 'F' && in_list(sel, "2F1"))
        return 1;

    return 0;
}

static int
class_selected(const char * tags)
{
    const char * s = opt.class_sel;

    if (s == NULL)
        return 1;

    while (1)
    {
        /* a term a+b+^c up to the next comma */
        int all = 1;

        while (*s != '\0' && *s != ',')
        {
            const char * e = s;
            int neg = 0, has;

            if (*s == '^')
            {
                neg = 1;
                s++;
                e = s;
            }

            while (*e != '\0' && *e != ',' && *e != '+')
                e++;

            has = has_tag(tags, s, e - s);
            if (has == neg)
                all = 0;

            s = e;
            if (*s == '+')
                s++;
        }

        if (all)
            return 1;
        if (*s == '\0')
            return 0;
        s++;
    }
}

static int
case_selected(const case_struct * cs)
{
    char id[2];

    if (!func_selected(cs->func) || !class_selected(cs->tags))
        return 0;

    if (opt.case_sel == NULL)
        return 1;

    id[0] = cs->func->id;
    id[1] = '\0';

    return in_list(opt.case_sel, cs->id) || in_list(opt.case_sel, id);
}

static void
print_call(const case_struct * cs)
{
    int i;

    printf("%s(", cs->func->name);
    for (i = 0; i < cs->func->nargs; i++)
        printf("%s%s = %s", i ? ", " : "", cs->func->argnames[i], cs->entry->args[i]);
    printf(")");
}

/* ------------------------------------------------------------------------- */
/*  Measurements                                                             */
/* ------------------------------------------------------------------------- */

#define STATUS_OK 0
#define STATUS_NA 1         /* the algorithm does not apply */
#define STATUS_PARSE 2
#define STATUS_TIMEOUT 3
#define STATUS_CRASH 4

typedef struct
{
    int status;
    int finite;
    int exact;
    int signal;
    int steps;              /* loops: the number of precisions */
    int reached;            /* loops: the goal was reached */
    int stalled;            /* loops: stopped as the inputs allow no more */
    int inputs_exact;       /* the inputs are exact at the (last) precision */
    slong acc;
    slong wp;               /* the last precision */
    double time;
    double first;           /* the time of the first call */
}
meas_struct;

typedef struct
{
    const case_struct * cs;
    const alg_struct * alg;
    int loop;
    slong prec;
    double timeout;
}
task_struct;

static slong
accuracy(const acb_t res, int * finite, int * exact)
{
    slong acc;

    *finite = acb_is_finite(res);
    *exact = 0;

    if (!*finite)
        return 0;

    if (acb_is_exact(res))
    {
        *exact = 1;
        return WORD_MAX;
    }

    acc = acb_rel_accuracy_bits(res);
    return FLINT_MAX(acc, 0);
}

static int
eval_inputs(acb_ptr x, const case_struct * cs, slong prec)
{
    int i, ok = 1;

    for (i = 0; i < cs->func->nargs; i++)
        if (!eval_expr(x + i, cs->entry->args[i], prec))
            ok = 0;

    return ok;
}

/* Whether some input is less accurate than the precision by more than
   16 bits, so that a higher precision need not improve the result */
static int
inputs_coarse(acb_srcptr x, int n, slong prec)
{
    int i;

    for (i = 0; i < n; i++)
        if (!acb_is_exact(x + i) && acb_rel_accuracy_bits(x + i) < prec - 16)
            return 1;

    return 0;
}

static void
measure(meas_struct * m, acb_t res, const task_struct * task)
{
    const case_struct * cs = task->cs;
    eval_func f = task->alg->f;
    int n = cs->func->nargs;
    acb_ptr x[MAX_LOOP_STEPS];
    slong wps[MAX_LOOP_STEPS];
    slong wp, prev_acc = 0, reps, r, goal;
    int k, steps = 0, prev_coarse = 0;
    double t, el, limit;
    acb_t tmp;

    memset(m, 0, sizeof(meas_struct));
    acb_indeterminate(res);
    acb_init(tmp);

    if (opt.threads > 1)
        flint_set_num_threads(opt.threads);

    goal = task->prec - opt.slack;
    m->first = 0.0;

    for (wp = task->prec; ; wp *= 2)
    {
        int ok;

        x[steps] = _acb_vec_init(n);
        wps[steps] = wp;
        steps++;

        if (!eval_inputs(x[steps - 1], cs, wp))
        {
            m->status = STATUS_PARSE;
            goto cleanup;
        }

        t = now();
        ok = f(res, x[steps - 1], wp);
        m->first += now() - t;

        if (!ok)
        {
            m->status = STATUS_NA;
            goto cleanup;
        }

        m->acc = accuracy(res, &m->finite, &m->exact);
        m->wp = wp;
        m->inputs_exact = 1;
        for (k = 0; k < n; k++)
            if (!acb_is_exact(x[steps - 1] + k))
                m->inputs_exact = 0;

        if (!task->loop)
            break;

        if (m->finite && m->acc >= goal)
        {
            m->reached = 1;
            break;
        }

        if (inputs_coarse(x[steps - 1], n, wp))
        {
            if (prev_coarse && m->acc <= prev_acc + 2)
            {
                m->stalled = 1;
                break;
            }
            prev_coarse = 1;
            prev_acc = m->acc;
        }

        if (2 * wp > opt.maxwp * task->prec || steps == MAX_LOOP_STEPS)
            break;
    }

    m->steps = steps;
    m->time = m->first;
    m->status = STATUS_OK;

    /* Repeat, unless the first call was slow: the least time of three
       batches which take a third of the minimum time each */
    limit = FLINT_MIN(0.5, task->timeout / 4);

    if (m->first < limit)
    {
        int batch;

        reps = 1;

        for (batch = 0; batch < 3; )
        {
            t = now();
            for (r = 0; r < reps; r++)
                for (k = 0; k < steps; k++)
                    f(tmp, x[k], wps[k]);
            el = now() - t;

            if (batch == 0 && 3 * el < opt.mintime && el < limit / 4)
            {
                reps *= (48 * el < opt.mintime) ? 8 : 2;
                continue;
            }

            if (batch == 0 || el / reps < m->time)
                m->time = el / reps;

            if (el >= limit / 4)
                break;

            batch++;
        }
    }

    else if (opt.warm && 3.0 * m->first < task->timeout)
    {
        /* once more, with the caches filled */
        t = now();
        for (k = 0; k < steps; k++)
            f(tmp, x[k], wps[k]);
        el = now() - t;

        if (el < m->time)
            m->time = el;
    }

cleanup:
    for (k = 0; k < steps; k++)
        _acb_vec_clear(x[k], n);
    acb_clear(tmp);
}

static void
work_measure(buf_t out, void * arg)
{
    meas_struct m;
    acb_t res;

    acb_init(res);
    measure(&m, res, (const task_struct *) arg);
    buf_put(out, &m, sizeof(m));
    buf_put_acb(out, res);
    acb_clear(res);
}

static void
run_measure(meas_struct * m, acb_t res, const task_struct * task)
{
    buf_t out;
    int status, sig;

    buf_init(out);
    status = run_isolated(out, work_measure, (void *) task, task->timeout, &sig);

    if (status == RUN_OK && buf_get(out, m, sizeof(meas_struct)) && buf_get_acb(out, res))
    {
        /* done */
    }
    else
    {
        memset(m, 0, sizeof(meas_struct));
        m->status = (status == RUN_TIMEOUT) ? STATUS_TIMEOUT : STATUS_CRASH;
        m->signal = sig;
        acb_indeterminate(res);
    }

    buf_clear(out);
}

/* Reference values: the function at the midpoints of the inputs and at
   the two points on the boundary in each real coordinate which has a
   radius, to prec + 24 bits, from the default function */

typedef struct
{
    int n;
    int ok[MAX_REF_POINTS];
    acb_struct v[MAX_REF_POINTS];
}
ref_struct;

static int
ref_value(acb_t v, const func_struct * func, acb_srcptr pt, slong prec)
{
    slong wp;
    int step, finite, exact;

    for (wp = prec + 64, step = 0; step < 7; step++, wp *= 2)
    {
        func->algs[0].f(v, pt, wp);
        if (accuracy(v, &finite, &exact) >= prec + 24 && finite)
            return 1;
    }

    return 0;
}

static void
work_ref(buf_t out, void * arg)
{
    const task_struct * task = arg;
    const case_struct * cs = task->cs;
    int n = cs->func->nargs, npts = 0, i, part, sign, ok;
    acb_ptr x, pt;
    acb_t v;
    arf_t r;

    x = _acb_vec_init(n);
    pt = _acb_vec_init(n);
    acb_init(v);
    arf_init(r);

    if (!eval_inputs(x, cs, task->prec))
    {
        buf_put(out, &npts, sizeof(int));
    }
    else
    {
        npts = 1;
        for (i = 0; i < n; i++)
            for (part = 0; part < 2; part++)
                if (!mag_is_zero(arb_radref(part ? acb_imagref(x + i) : acb_realref(x + i))))
                    npts += 2;

        buf_put(out, &npts, sizeof(int));

        for (i = 0; i < n; i++)
            acb_get_mid(pt + i, x + i);

        ok = ref_value(v, cs->func, pt, task->prec);
        buf_put(out, &ok, sizeof(int));
        buf_put_acb(out, v);

        for (i = 0; i < n && ok; i++)
        {
            for (part = 0; part < 2; part++)
            {
                arb_srcptr c = part ? acb_imagref(x + i) : acb_realref(x + i);
                arb_ptr p = part ? acb_imagref(pt + i) : acb_realref(pt + i);
                int ok2;

                if (mag_is_zero(arb_radref(c)))
                    continue;

                for (sign = 0; sign < 2; sign++)
                {
                    arf_set_mag(r, arb_radref(c));
                    if (sign)
                        arf_neg(r, r);
                    arf_add(arb_midref(p), arb_midref(c), r, ARF_PREC_EXACT, ARF_RND_DOWN);
                    ok2 = ref_value(v, cs->func, pt, task->prec);
                    buf_put(out, &ok2, sizeof(int));
                    buf_put_acb(out, v);
                }

                arf_set(arb_midref(p), arb_midref(c));
            }
        }
    }

    _acb_vec_clear(x, n);
    _acb_vec_clear(pt, n);
    acb_clear(v);
    arf_clear(r);
}

/* Returns the accuracy which the inputs allow (in bits; HUGE_VAL if the
   inputs are exact, and a negative number if it is not known). */
static double
run_ref(ref_struct * ref, const case_struct * cs, slong prec)
{
    task_struct task;
    buf_t out;
    int status, sig, i, npts;
    double attain = -1.0;

    task.cs = cs;
    task.alg = cs->func->algs;
    task.loop = 0;
    task.prec = prec;
    task.timeout = 4 * opt.timeout;

    ref->n = 0;
    buf_init(out);
    status = run_isolated(out, work_ref, &task, task.timeout, &sig);

    if (status == RUN_OK && buf_get(out, &npts, sizeof(int)) && npts >= 1 && npts <= MAX_REF_POINTS)
    {
        for (i = 0; i < npts; i++)
        {
            if (!buf_get(out, &ref->ok[i], sizeof(int)) || !buf_get_acb(out, ref->v + i))
                break;
            ref->n = i + 1;
        }

        if (ref->n >= 1 && !ref->ok[0])
            ref->n = 0;

        if (ref->n == npts)
        {
            /* |f(mid)| over the sum over the coordinates of the larger
               of the two differences */
            mag_t sum, d, e;
            acb_t t;
            int all = 1;

            mag_init(sum);
            mag_init(d);
            mag_init(e);
            acb_init(t);

            for (i = 1; i + 1 < npts; i += 2)
            {
                if (!ref->ok[i] || !ref->ok[i + 1])
                {
                    all = 0;
                    break;
                }
                acb_sub(t, ref->v + i, ref->v, prec + 64);
                acb_get_mag(d, t);
                acb_sub(t, ref->v + i + 1, ref->v, prec + 64);
                acb_get_mag(e, t);
                mag_max(d, d, e);
                mag_add(sum, sum, d);
            }

            if (all)
            {
                if (npts == 1)
                    attain = HUGE_VAL;
                else
                {
                    acb_get_mag_lower(d, ref->v);
                    if (mag_is_zero(d))
                        attain = -1.0;
                    else if (mag_is_zero(sum))
                        attain = (double) prec;
                    else
                    {
                        attain = mag_get_d_log2_approx(d) - mag_get_d_log2_approx(sum);
                        attain = FLINT_MAX(attain, 0.0);
                        attain = FLINT_MIN(attain, (double) prec);
                    }
                }
            }

            mag_clear(sum);
            mag_clear(d);
            mag_clear(e);
            acb_clear(t);
        }
    }

    buf_clear(out);
    return attain;
}

/* ------------------------------------------------------------------------- */
/*  Records and output                                                       */
/* ------------------------------------------------------------------------- */

#define FLAG_FASTER 1           /* as accurate as the default and faster */
#define FLAG_ACCURATE 2         /* more accurate than the default */
#define FLAG_INCONSISTENT 4     /* does not overlap another result */
#define FLAG_CONTAIN 8          /* does not contain a reference value */
#define FLAG_LESS 16            /* the default: less accurate than the goal */
#define FLAG_UNJUDGED 32        /* the default: less accurate, inexact inputs, no -ref */

typedef struct
{
    slong cs;
    int prec_idx;
    int alg;
    int loop;
    int flags;
    meas_struct m;
    double attain;
    double ratio;       /* time over the time of the default; 0 if not known */
    slong def_acc;      /* the accuracy of the default (-1 if no result) */
    double def_time;
    double limit;       /* the time limit */
}
record_struct;

static record_struct * records;
static slong num_records, alloc_records;

/* A case and precision at which some algorithm or loop reaches the goal
   (the accuracy which counts as full accuracy) faster than the default
   function does: directly or, if it is less accurate, in its loop */
typedef struct
{
    slong def_rec;      /* the record of the default, or -1 if it does not reach the goal */
    slong alt_rec;      /* the fastest other record which reaches the goal */
    slong first_rec;    /* the record of the direct call of the default */
    int loop_tried;     /* the loop of the default was run */
    double ratio;       /* time over the time of the default; 0 if def_rec is -1 */
}
route_struct;

static route_struct * routes;
static slong num_routes, alloc_routes;

static int
reaches_goal(const record_struct * rec, double goal)
{
    return rec->m.status == STATUS_OK && rec->m.finite && goal > 0
        && (rec->m.exact || rec->m.acc >= goal);
}

/* Compares the records first, ..., num_records - 1 of one case and precision */
static void
find_route(slong first, slong prec)
{
    double goal;
    slong i, def = -1, alt = -1;
    int loop_tried = 0;

    goal = records[first].attain;
    goal = (goal >= 0.0 && goal < prec) ? goal : (double) prec;
    goal -= opt.slack;

    if (reaches_goal(records + first, goal))
        def = first;

    for (i = first + 1; i < num_records; i++)
    {
        if (records[i].alg == 0)
        {
            loop_tried = 1;
            if (def < 0 && reaches_goal(records + i, goal))
                def = i;
        }
        else if (reaches_goal(records + i, goal)
            && (alt < 0 || records[i].m.time < records[alt].m.time))
        {
            alt = i;
        }
    }

    if (alt >= 0 && (def < 0 || records[alt].m.time * opt.ratio < records[def].m.time))
    {
        route_struct * r;

        if (num_routes == alloc_routes)
        {
            alloc_routes = 2 * alloc_routes + 64;
            routes = flint_realloc(routes, alloc_routes * sizeof(route_struct));
        }

        r = routes + num_routes++;
        r->def_rec = def;
        r->alt_rec = alt;
        r->first_rec = first;
        r->loop_tried = loop_tried;
        r->ratio = (def >= 0 && records[def].m.time > 0.0) ?
            records[alt].m.time / records[def].m.time : 0.0;
    }
}

static record_struct *
new_record(void)
{
    if (num_records == alloc_records)
    {
        alloc_records = 2 * alloc_records + 256;
        records = flint_realloc(records, alloc_records * sizeof(record_struct));
    }

    memset(records + num_records, 0, sizeof(record_struct));
    return records + num_records++;
}

static void
format_time(char * s, size_t n, double t)
{
    if (t < 999.5e-9)
        snprintf(s, n, "%.3g ns", t * 1e9);
    else if (t < 999.5e-6)
        snprintf(s, n, "%.3g us", t * 1e6);
    else if (t < 999.5e-3)
        snprintf(s, n, "%.3g ms", t * 1e3);
    else
        snprintf(s, n, "%.3g s", t);
}

static void
alg_label(char * s, size_t n, const func_struct * func, int alg, int loop)
{
    snprintf(s, n, "%s%s", loop ? "loop:" : "", func->algs[alg].name);
}

static int
usable(const meas_struct * m)
{
    return m->status == STATUS_OK && m->finite;
}

/* The accuracy which counts as full accuracy */
static double
goal_bits(slong prec, double attain)
{
    double g = (double) prec;
    if (attain >= 0.0 && attain < g)
        g = attain;
    return g - opt.slack;
}

static void
print_value(const char * label, const acb_t x, slong prec)
{
    slong digits = FLINT_MIN(prec, 128) * 0.30103 + 1;
    printf("                 %s = ", label);
    acb_printn(x, digits, 0);
    printf("\n");
}

static void
csv_string(FILE * f, const char * s)
{
    fputc('"', f);
    for ( ; *s != '\0'; s++)
    {
        if (*s == '"')
            fputc('"', f);
        fputc(*s, f);
    }
    fputc('"', f);
}

static const char *
status_name(int status)
{
    switch (status)
    {
        case STATUS_OK: return "ok";
        case STATUS_NA: return "na";
        case STATUS_PARSE: return "invalid";
        case STATUS_TIMEOUT: return "timeout";
        default: return "aborted";
    }
}

static void
csv_header(FILE * f)
{
    fprintf(f, "case,function,inputs,classes,prec,algorithm,status,time,first_call,accuracy,"
        "exact,last_prec,steps,reached,attainable,time_ratio,flags\n");
}

static void
csv_record(FILE * f, const record_struct * rec)
{
    const case_struct * cs = cases + rec->cs;
    const meas_struct * m = &rec->m;
    char label[64], inputs[1024];
    size_t len = 0;
    int i;

    alg_label(label, sizeof(label), cs->func, rec->alg, rec->loop);

    inputs[0] = '\0';
    for (i = 0; i < cs->func->nargs && len < sizeof(inputs); i++)
        len += snprintf(inputs + len, sizeof(inputs) - len, "%s%s = %s", i ? "; " : "",
            cs->func->argnames[i], cs->entry->args[i]);

    fprintf(f, "%s,%s,", cs->id, cs->func->name);
    csv_string(f, inputs);
    fputc(',', f);
    csv_string(f, cs->tags);
    fprintf(f, ",%ld,%s,%s,", (long) opt.precs[rec->prec_idx], label, status_name(m->status));

    if (m->status == STATUS_OK)
    {
        fprintf(f, "%.6g,%.6g,", m->time, m->first);
        if (m->finite)
            fprintf(f, "%ld", (long) (m->exact ? m->wp : m->acc));
        fprintf(f, ",%d,%ld,%d,%s,", m->exact, (long) m->wp, m->steps,
            rec->loop ? (m->reached ? "1" : "0") : "");
    }
    else
    {
        fprintf(f, ",,,,,,,");
    }

    if (rec->attain == HUGE_VAL)
        fprintf(f, "inf");
    else if (rec->attain >= 0.0)
        fprintf(f, "%.0f", rec->attain);
    fputc(',', f);

    if (rec->ratio > 0.0)
        fprintf(f, "%.4g", rec->ratio);

    fprintf(f, ",%s%s%s%s%s\n",
        (rec->flags & FLAG_LESS) ? "less-accurate " : "",
        (rec->flags & FLAG_FASTER) ? "faster " : "",
        (rec->flags & FLAG_ACCURATE) ? "more-accurate " : "",
        (rec->flags & FLAG_INCONSISTENT) ? "INCONSISTENT " : "",
        (rec->flags & FLAG_CONTAIN) ? "NOT-CONTAINED " : "");
}

static void
print_record(const record_struct * rec, const acb_t res)
{
    const case_struct * cs = cases + rec->cs;
    const meas_struct * m = &rec->m;
    slong prec = opt.precs[rec->prec_idx];
    char label[64], tstr[32], accstr[32], attstr[32];

    alg_label(label, sizeof(label), cs->func, rec->alg, rec->loop);

    tstr[0] = accstr[0] = attstr[0] = '\0';

    if (m->status == STATUS_OK)
    {
        format_time(tstr, sizeof(tstr), m->time);

        if (!m->finite)
            snprintf(accstr, sizeof(accstr), "none");
        else if (m->exact)
            snprintf(accstr, sizeof(accstr), "exact");
        else
            snprintf(accstr, sizeof(accstr), "%ld", (long) m->acc);
    }

    if (rec->attain == HUGE_VAL)
        snprintf(attstr, sizeof(attstr), "-");
    else if (rec->attain >= 0.0)
        snprintf(attstr, sizeof(attstr), "%.0f", rec->attain);
    else
        snprintf(attstr, sizeof(attstr), "?");

    printf("    %7ld  %-20s %10s %7s", (long) prec, label, tstr, accstr);
    if (opt.ref)
        printf(" %7s", (rec->alg == 0 && !rec->loop) ? attstr : "");

    if (m->status == STATUS_NA)
        printf("  not applicable");
    else if (m->status == STATUS_TIMEOUT)
        printf("  stopped after %.3g s", rec->limit);
    else if (m->status == STATUS_CRASH)
        printf("  aborted (signal %d)", m->signal);
    else if (m->status == STATUS_PARSE)
        printf("  invalid input");

    if (m->status == STATUS_OK && rec->loop)
    {
        printf("  wp %ld", (long) m->wp);
        if (!m->reached)
            printf(m->stalled ? " (inputs allow no more)" : " (goal not reached)");
    }

    if (m->status == STATUS_OK && m->first > 1e-3 && m->first > 5.0 * m->time)
    {
        char fstr[32];
        format_time(fstr, sizeof(fstr), m->first);
        printf("  (first call %s)", fstr);
    }

    if (rec->ratio > 0.0)
        printf("  %.3gx", rec->ratio);
    if (rec->flags & FLAG_FASTER)
        printf("  <-- faster");
    if (rec->flags & FLAG_ACCURATE)
        printf("  <-- more accurate");
    if (rec->flags & FLAG_INCONSISTENT)
        printf("  INCONSISTENT");
    if (rec->flags & FLAG_CONTAIN)
        printf("  DOES NOT CONTAIN A REFERENCE VALUE");
    printf("\n");

    if (opt.values && m->status == STATUS_OK && res != NULL)
        print_value("", res, rec->loop ? m->wp : prec);
}

/* The algorithms to run for a function: plan[k] = alg + MAX_ALGS * loop */
static int
build_plan(int * plan, const func_struct * func)
{
    int n = 0, a, loop;
    const char * sel = opt.alg_sel;
    char name[80];

    plan[n++] = 0;

    if (sel == NULL)
        return n;

    for (loop = 0; loop < 2; loop++)
    {
        for (a = 0; func->algs[a].name != NULL && a < MAX_ALGS; a++)
        {
            int want;

            if (!loop)
            {
                if (a == 0)
                    continue;
                want = in_list(sel, "all") || in_list(sel, func->algs[a].name);
            }
            else
            {
                snprintf(name, sizeof(name), "loop:%s", func->algs[a].name);
                want = in_list(sel, "loop:all") || in_list(sel, name)
                    || (a == 0 && in_list(sel, "loop"));
            }

            if (want)
                plan[n++] = a + MAX_ALGS * loop;
        }
    }

    return n;
}

/* Checks that the algorithm names of -alg exist */
static int
check_alg_names(void)
{
    const char * s = opt.alg_sel;
    int ok = 1;

    while (s != NULL && *s != '\0')
    {
        const char * e = strchr(s, ',');
        size_t n = (e == NULL) ? strlen(s) : (size_t) (e - s);
        char name[80];
        const char * p = name;
        int found = 0;
        slong i, a;

        if (n >= sizeof(name))
            n = sizeof(name) - 1;
        memcpy(name, s, n);
        name[n] = '\0';

        if (strncmp(name, "loop:", 5) == 0)
            p = name + 5;

        if (strcmp(p, "all") == 0 || strcmp(name, "loop") == 0)
            found = 1;

        for (i = 0; i < NUM_FUNCS && !found; i++)
            for (a = 0; funcs[i].algs[a].name != NULL; a++)
                if (strcmp(funcs[i].algs[a].name, p) == 0)
                    found = 1;

        if (!found)
        {
            printf("error: unknown algorithm \"%s\" (see -algs)\n", name);
            ok = 0;
        }

        if (e == NULL)
            break;
        s = e + 1;
    }

    return ok;
}

/* Measures the default function and the selected algorithms and loops
   for one case, checks that all results overlap, and prints them. */
static void
run_case(slong ci)
{
    const case_struct * cs = cases + ci;
    const func_struct * func = cs->func;
    int plan[2 * MAX_ALGS + 1], source[2 * MAX_ALGS + 1], nplan, p, k, i;
    acb_ptr results, block, x;
    slong num_results = 0, j;
    ref_struct ref;
    acb_t res;

    nplan = build_plan(plan, func);
    results = _acb_vec_init(nplan * opt.nprecs);
    block = _acb_vec_init(nplan);
    x = _acb_vec_init(func->nargs);
    acb_init(res);
    for (i = 0; i < MAX_REF_POINTS; i++)
        acb_init(ref.v + i);

    if (!opt.summary_only)
    {
        printf("%s  ", cs->id);
        print_call(cs);
        printf("\n      [%s]\n", cs->tags);
    }

    for (p = 0; p < opt.nprecs; p++)
    {
        slong prec = opt.precs[p];
        slong first = num_records;
        record_struct * recs;
        double attain = -1.0, goal;
        task_struct task;

        ref.n = 0;

        if (opt.values && !opt.summary_only && eval_inputs(x, cs, prec))
            for (i = 0; i < func->nargs; i++)
                print_value(func->argnames[i], x + i, prec);

        if (opt.ref)
            attain = run_ref(&ref, cs, prec);

        goal = goal_bits(prec, attain);

        for (k = 0; k < nplan; k++)
            new_record();

        recs = records + first;
        task.cs = cs;
        task.prec = prec;

        /* The measurements */
        for (k = 0; k < nplan; k++)
        {
            record_struct * rec = recs + k;

            rec->cs = ci;
            rec->prec_idx = p;
            rec->alg = plan[k] % MAX_ALGS;
            rec->loop = plan[k] / MAX_ALGS;
            rec->attain = attain;
            source[k] = -1;

            task.alg = func->algs + rec->alg;
            task.loop = rec->loop;
            task.timeout = opt.timeout;
            if (k > 0 && recs[0].m.status == STATUS_OK)
                task.timeout = FLINT_MIN(opt.timeout, FLINT_MAX(0.1, 100 * recs[0].m.first));
            rec->limit = task.timeout;

            /* The loop of an algorithm whose direct call is accurate is
               that call */
            for (i = 0; i < k && rec->loop; i++)
                if (recs[i].alg == rec->alg && !recs[i].loop && recs[i].m.status == STATUS_OK
                    && recs[i].m.finite && (recs[i].m.exact || recs[i].m.acc >= prec - opt.slack))
                    source[k] = i;

            if (source[k] < 0)
                run_measure(&rec->m, block + k, &task);
        }

        /* The default function once more, after the others: the lesser
           of its two times is compared with theirs */
        if (nplan > 1 && recs[0].m.status == STATUS_OK)
        {
            meas_struct again;

            task.alg = func->algs;
            task.loop = 0;
            task.timeout = opt.timeout;
            run_measure(&again, res, &task);

            if (again.status == STATUS_OK && again.time < recs[0].m.time)
                recs[0].m.time = again.time;
        }

        /* Comparisons, checks and output */
        for (k = 0; k < nplan; k++)
        {
            record_struct * rec = recs + k;
            meas_struct * m = &rec->m;
            const meas_struct * def = &recs[0].m;

            if (source[k] >= 0)
            {
                *m = recs[source[k]].m;
                m->reached = 1;
            }

            rec->def_acc = usable(def) ? (def->exact ? prec : def->acc) : -1;
            rec->def_time = def->time;

            if (k == 0)
            {
                /* Less accurate than the precision, or than the inputs
                   allow: this is known for exact inputs and with -ref */
                if (m->status == STATUS_OK && !m->exact && goal > 0
                    && (!m->finite || m->acc < goal))
                    rec->flags |= (m->inputs_exact || attain >= 0.0) ? FLAG_LESS : FLAG_UNJUDGED;
            }
            else if (usable(m) && def->status == STATUS_OK)
            {
                slong acc = m->exact ? prec : m->acc;

                rec->ratio = (def->time > 0.0) ? m->time / def->time : 0.0;

                if (rec->def_acc > 0 && acc >= FLINT_MIN((double) rec->def_acc, goal)
                    && m->time * opt.ratio < def->time)
                    rec->flags |= FLAG_FASTER;

                if (acc >= rec->def_acc + opt.slack && rec->def_acc < goal)
                    rec->flags |= FLAG_ACCURATE;
            }

            if (usable(m) && source[k] < 0)
            {
                for (j = 0; j < num_results; j++)
                    if (!acb_overlaps(block + k, results + j))
                        rec->flags |= FLAG_INCONSISTENT;

                /* (the results of the loops are for other input balls) */
                for (i = 0; i < ref.n && !rec->loop; i++)
                    if (ref.ok[i] && !acb_overlaps(block + k, ref.v + i))
                        rec->flags |= FLAG_CONTAIN;

                acb_set(results + num_results, block + k);
                num_results++;
            }

            if (!opt.summary_only)
                print_record(rec, source[k] < 0 ? block + k : NULL);

            if (opt.csv != NULL)
                csv_record(opt.csv, rec);
        }

        if (opt.csv != NULL)
            fflush(opt.csv);
        fflush(stdout);

        if (nplan > 1)
            find_route(first, prec);
    }

    if (!opt.summary_only)
        printf("\n");

    _acb_vec_clear(results, nplan * opt.nprecs);
    _acb_vec_clear(block, nplan);
    _acb_vec_clear(x, func->nargs);
    acb_clear(res);
    for (i = 0; i < MAX_REF_POINTS; i++)
        acb_clear(ref.v + i);
}

/* ------------------------------------------------------------------------- */
/*  Summary                                                                  */
/* ------------------------------------------------------------------------- */

/* The records with a flag */
static void
print_flagged(const char * title, int flag)
{
    slong n = 0, i;
    char label[64], tstr[32];

    for (i = 0; i < num_records; i++)
        n += (records[i].flags & flag) != 0;

    if (n == 0)
        return;

    printf("\n%s (%ld):\n\n", title, (long) n);

    for (i = 0; i < num_records; i++)
    {
        const record_struct * rec = records + i;
        const case_struct * cs = cases + rec->cs;

        if (!(rec->flags & flag))
            continue;

        printf("  %-5s %6ld  ", cs->id, (long) opt.precs[rec->prec_idx]);

        if (flag & FLAG_LESS)
        {
            format_time(tstr, sizeof(tstr), rec->m.time);
            printf("%9s ", tstr);
            if (rec->m.finite)
                printf("%6ld bits", (long) rec->m.acc);
            else
                printf("  no result");
            if (rec->attain >= 0.0 && rec->attain != HUGE_VAL)
                printf(" of %-6.0f ", rec->attain);
            else
                printf("           ");
        }
        else
        {
            alg_label(label, sizeof(label), cs->func, rec->alg, rec->loop);
            printf("%-20s ", label);
        }

        print_call(cs);
        printf("\n");
    }
}

static int
cmp_route(const void * a, const void * b)
{
    const route_struct * x = a;
    const route_struct * y = b;
    if (x->ratio < y->ratio)
        return -1;
    if (x->ratio > y->ratio)
        return 1;
    return (x->first_rec < y->first_rec) ? -1 : (x->first_rec > y->first_rec);
}

/* The cases in which the default function is not the fastest way to an
   accurate result */
static void
print_routes(void)
{
    slong i;
    char label[64], t1[32], t2[32], acc[32];

    if (num_routes == 0)
        return;

    qsort(routes, num_routes, sizeof(route_struct), cmp_route);

    printf("\nFaster ways to an accurate result than the default function (%ld):\n\n"
        "  The default function is taken alone if its result is accurate, and in its\n"
        "  loop otherwise (if the loop is run: -alg loop); it is compared with the\n"
        "  fastest other algorithm or loop with an accurate result, and listed if\n"
        "  that takes less than 1/%.3g of the time (time, and accuracy in bits).\n\n",
        (long) num_routes, opt.ratio);

    for (i = 0; i < num_routes; i++)
    {
        const route_struct * r = routes + i;
        const record_struct * alt = records + r->alt_rec;
        const record_struct * def = records + (r->def_rec >= 0 ? r->def_rec : r->first_rec);
        const case_struct * cs = cases + alt->cs;
        slong prec = opt.precs[alt->prec_idx];

        format_time(t1, sizeof(t1), def->m.time);
        format_time(t2, sizeof(t2), alt->m.time);

        if (def->m.status != STATUS_OK)
            snprintf(acc, sizeof(acc), "failed");
        else if (!def->m.finite)
            snprintf(acc, sizeof(acc), "none");
        else
            snprintf(acc, sizeof(acc), "%ld", (long) (def->m.exact ? prec : def->m.acc));

        printf("  %-5s %6ld  ", cs->id, (long) prec);

        if (r->def_rec >= 0)
        {
            alg_label(label, sizeof(label), cs->func, def->alg, def->loop);
            printf("%-13s %9s %6s   ", label, t1, acc);
        }
        else
        {
            printf("%-13s %9s %6s   ", r->loop_tried ? "not in loop" : "not accurate",
                def->m.status == STATUS_OK ? t1 : "", acc);
        }

        alg_label(label, sizeof(label), cs->func, alt->alg, alt->loop);
        printf("%-20s %9s %6ld   ", label, t2, (long) (alt->m.exact ? prec : alt->m.acc));

        if (r->def_rec >= 0)
            printf("%6.3gx   ", r->ratio);
        else
            printf("          ");

        print_call(cs);
        printf("\n");
    }
}

static void
print_summary(double elapsed)
{
    slong f, i;
    int p;
    char tstr[32];

    printf("Summary of the default functions\n\n");
    if (opt.ref)
        printf("  accurate: within %ld bits of the precision or of what the inputs allow\n"
               "  less:     finite and less accurate;  none: not finite;  failed: stopped\n\n",
               (long) opt.slack);
    else
        printf("  accurate: within %ld bits of the precision\n"
               "  less:     less accurate or not finite, with exact inputs\n"
               "  inexact:  less accurate or not finite, with inexact inputs (use -ref\n"
               "            to compare with the accuracy which the inputs allow)\n"
               "  failed:   stopped\n\n", (long) opt.slack);
    printf("  %-16s %6s %6s %9s %7s %7s %7s %11s", "function", "prec", "cases",
        "accurate", "less", opt.ref ? "none" : "inexact", "failed", "total time");
    if (opt.alg_sel != NULL)
        printf(" %8s %9s", "slower", "only alg.");
    printf("\n");

    for (f = 0; f < NUM_FUNCS; f++)
    {
        for (p = 0; p < opt.nprecs; p++)
        {
            slong n = 0, good = 0, less = 0, none = 0, failed = 0, slower = 0, noroute = 0;
            double total = 0.0;

            for (i = 0; i < num_records; i++)
            {
                const record_struct * rec = records + i;

                if (cases[rec->cs].func != funcs + f || rec->prec_idx != p)
                    continue;

                if (rec->alg == 0 && !rec->loop)
                {
                    n++;
                    if (rec->m.status != STATUS_OK)
                        failed++;
                    else if (rec->flags & FLAG_UNJUDGED)
                        none++;
                    else if (!(rec->flags & FLAG_LESS))
                        good++;
                    else if (opt.ref && !rec->m.finite)
                        none++;
                    else
                        less++;

                    if (rec->m.status == STATUS_OK)
                        total += rec->m.time;
                }
            }

            for (i = 0; i < num_routes; i++)
            {
                const record_struct * rec = records + routes[i].first_rec;

                if (cases[rec->cs].func == funcs + f && rec->prec_idx == p)
                {
                    if (routes[i].def_rec >= 0)
                        slower++;
                    else
                        noroute++;
                }
            }

            if (n == 0)
                continue;

            format_time(tstr, sizeof(tstr), total);
            printf("  %-16s %6ld %6ld %9ld %7ld %7ld %7ld %11s", funcs[f].name,
                (long) opt.precs[p], (long) n, (long) good, (long) less, (long) none,
                (long) failed, tstr);
            if (opt.alg_sel != NULL)
                printf(" %8ld %9ld", (long) slower, (long) noroute);
            printf("\n");
        }
    }

    print_flagged(opt.ref ?
        "The default function is less accurate than the precision and the inputs allow" :
        "The default function is less accurate than the precision, with exact inputs",
        FLAG_LESS);

    if (opt.alg_sel != NULL)
        print_routes();

    print_flagged("INCONSISTENT RESULTS (not overlapping another result of the same case)",
        FLAG_INCONSISTENT);
    print_flagged("RESULTS WHICH DO NOT CONTAIN A REFERENCE VALUE", FLAG_CONTAIN);

    printf("\n%ld measurements in %.1f s\n", (long) num_records, elapsed);
}

/* ------------------------------------------------------------------------- */
/*  Main                                                                     */
/* ------------------------------------------------------------------------- */

static int
parse_precs(const char * s)
{
    opt.nprecs = 0;

    while (*s != '\0')
    {
        char * end;
        long v = strtol(s, &end, 10);

        if (end == s || v < 2 || opt.nprecs == MAX_PRECS)
            return 0;
        opt.precs[opt.nprecs++] = v;
        s = end;
        if (*s == ',')
            s++;
        else if (*s != '\0')
            return 0;
    }

    return opt.nprecs > 0;
}

int
main(int argc, char * argv[])
{
    slong i, selected = 0;
    int a;
    double t0;

    opt.precs[0] = 64;
    opt.precs[1] = 333;
    opt.precs[2] = 3333;
    opt.nprecs = 3;
    opt.func_sel = opt.class_sel = opt.case_sel = opt.alg_sel = NULL;
    opt.ninputs = 0;
    opt.list = opt.list_algs = opt.ref = opt.values = 0;
    opt.csv_name = NULL;
    opt.csv = NULL;
    opt.summary_only = opt.verbose = 0;
    opt.use_fork = CATALOG_HAVE_FORK;
    opt.threads = 1;
    opt.mintime = 0.01;
    opt.timeout = 10.0;
    opt.ratio = 1.5;
    opt.slack = 10;
    opt.maxmem = 4096;
    opt.maxwp = 32;

    for (a = 1; a < argc; a++)
    {
        const char * o = argv[a];
        const char * v = (a + 1 < argc) ? argv[a + 1] : NULL;

        if (o[0] == '-' && o[1] == '-')
            o++;

        if (!strcmp(o, "-help") || !strcmp(o, "-h"))
        {
            usage();
            return 0;
        }
        else if (!strcmp(o, "-list"))
            opt.list = 1;
        else if (!strcmp(o, "-algs"))
            opt.list_algs = 1;
        else if (!strcmp(o, "-ref"))
            opt.ref = 1;
        else if (!strcmp(o, "-values"))
            opt.values = 1;
        else if (!strcmp(o, "-summary"))
            opt.summary_only = 1;
        else if (!strcmp(o, "-nofork"))
            opt.use_fork = 0;
        else if (!strcmp(o, "-verbose"))
            opt.verbose = 1;
        else if (!strcmp(o, "-warm"))
            opt.warm = 1;
        else if (v == NULL)
        {
            printf("error: unknown option or missing value: %s (see -help)\n", argv[a]);
            return 1;
        }
        else
        {
            a++;

            if (!strcmp(o, "-func"))
                opt.func_sel = v;
            else if (!strcmp(o, "-class"))
                opt.class_sel = v;
            else if (!strcmp(o, "-case"))
                opt.case_sel = v;
            else if (!strcmp(o, "-input"))
            {
                if (opt.ninputs == MAX_INPUTS)
                {
                    printf("error: too many inputs\n");
                    return 1;
                }
                opt.inputs[opt.ninputs++] = v;
            }
            else if (!strcmp(o, "-alg"))
                opt.alg_sel = v;
            else if (!strcmp(o, "-csv"))
                opt.csv_name = v;
            else if (!strcmp(o, "-prec"))
            {
                if (!parse_precs(v))
                {
                    printf("error: invalid list of precisions: %s\n", v);
                    return 1;
                }
            }
            else if (!strcmp(o, "-mintime"))
                opt.mintime = atof(v);
            else if (!strcmp(o, "-timeout"))
                opt.timeout = atof(v);
            else if (!strcmp(o, "-ratio"))
                opt.ratio = atof(v);
            else if (!strcmp(o, "-slack"))
                opt.slack = atol(v);
            else if (!strcmp(o, "-maxmem"))
                opt.maxmem = atol(v);
            else if (!strcmp(o, "-maxwp"))
                opt.maxwp = atol(v);
            else if (!strcmp(o, "-threads"))
                opt.threads = atoi(v);
            else
            {
                printf("error: unknown option: %s (see -help)\n", argv[a - 1]);
                return 1;
            }
        }
    }

    if (opt.timeout <= 0.0 || opt.mintime < 0.0 || opt.ratio < 1.0 || opt.maxwp < 1)
    {
        printf("error: invalid value of -timeout, -mintime, -ratio or -maxwp\n");
        return 1;
    }

    if (opt.list_algs)
    {
        for (i = 0; i < NUM_FUNCS; i++)
        {
            printf("%-16s %c   %s\n", funcs[i].name, funcs[i].id, funcs[i].desc);
            for (a = 0; funcs[i].algs[a].name != NULL; a++)
                printf("    %-16s %s\n", funcs[i].algs[a].name, funcs[i].algs[a].desc);
        }
        return 0;
    }

    if (!check_alg_names())
        return 1;

    if (!(opt.ninputs > 0 ? build_inputs() : build_cases()))
        return 1;

    for (i = 0; i < num_cases; i++)
        selected += case_selected(cases + i);

    if (selected == 0)
    {
        printf("error: no case is selected (see -list)\n");
        return 1;
    }

    if (opt.list)
    {
        for (i = 0; i < num_cases; i++)
        {
            if (case_selected(cases + i))
            {
                printf("%-5s ", cases[i].id);
                print_call(cases + i);
                printf("\n      [%s]\n", cases[i].tags);
            }
        }
        printf("\n%ld cases\n", (long) selected);
        return 0;
    }

    if (opt.csv_name != NULL)
    {
        opt.csv = fopen(opt.csv_name, "w");
        if (opt.csv == NULL)
        {
            printf("error: cannot write %s\n", opt.csv_name);
            return 1;
        }
        csv_header(opt.csv);
    }

    if (!opt.summary_only)
    {
        printf("       prec  algorithm                  time     acc");
        if (opt.ref)
            printf("  inputs");
        printf("\n\n");
    }

    t0 = now();

    for (i = 0; i < num_cases; i++)
        if (case_selected(cases + i))
            run_case(i);

    print_summary(now() - t0);

    if (opt.csv != NULL)
        fclose(opt.csv);

    flint_free(cases);
    flint_free(records);
    flint_free(routes);
    flint_cleanup_master();
    return 0;
}
