/* 
    Examples for double exponential integration demonstrating integration 
    with log, log^2 and sqrt singularities. In all cases close to
    quadratic convergence with the number of function evaluations is 
    observed. 

    This file is part of FLINT.

    Copyright (C) 2024 Hartmut Monien
*/

#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include <flint/profiler.h>
#include <flint/arb_hypgeom.h>
#include <flint/acb_hypgeom.h>
#include <flint/acb_dirichlet.h>
#include <flint/acb_modular.h>
#include <flint/arb_calc.h>

enum QUAD_RULE
{
    QUAD_INTERVAL,
    QUAD_SEMI_INF,
    QUAD_INF
};

void
quad_rule(arb_t w, arb_t t, arb_t x, arb_t pi, int rule, slong prec)
{

    arb_t c, s;

    arb_init(s);
    arb_init(c);

    arb_sinh_cosh(s, c, x, prec);
    arb_mul(c, c, pi, prec);
    arb_mul(s, s, pi, prec);
    arb_div_ui(s, s, 2, prec);
    arb_div_ui(c, c, 2, prec);

    switch (rule)
    {
        case QUAD_INTERVAL:
            arb_sinh_cosh(t, w, s, prec);
            arb_div(t, t, w, prec);
            arb_sqr(w, w, prec);
            arb_div(w, c, w, prec);
            break;
        case QUAD_SEMI_INF:
            arb_exp(t, s, prec);
            arb_mul(w, c, t, prec);
            break;
        case QUAD_INF:
            arb_set(w, c);
            arb_sinh_cosh(t, c, s, prec);
            arb_mul(w, w, c, prec);
            break;
        default:
            flint_printf("no integration rule %d.\n", rule);
            exit(1);
    }

    arb_clear(s);
    arb_clear(c);
}

void
quad(arb_t result,
     arb_calc_func_t f, void *param, int rule, int verbose, slong prec)
{

    arb_t x, y, h, x0, pi, w, t;
    slong order = 0, max_iter = 16, n_eval = 0, flag, j, k;

    arb_ptr r = _arb_vec_init(max_iter);
    arb_ptr d = _arb_vec_init(max_iter);

    arb_init(x);
    arb_init(y);
    arb_init(h);
    arb_init(x0);

    arb_init(pi);

    arb_init(w);
    arb_init(t);

    arb_const_pi(pi, prec);

    arb_zero(t);
    quad_rule(w, t, x, pi, rule, prec);

    flag = f(y, t, param, order, prec);
    n_eval++;
    arb_mul(result, y, w, prec);
    arb_set(r, result);

    if (verbose > 0)
    {
        flint_printf("\n\t%8s %8s %10s\n\n", "log(d)", "#evals", "rate");
    }

    arb_one(x0);
    arb_one(h);

    for (j = 0; j + 1 < max_iter; j++)
    {

        for (int pm = -1; pm <= 1; pm += 2)
        {

            arb_set(x, x0);

            if (pm < 0)
            {
                arb_neg(x, x);
            }

            for (k = 0;; k++)
            {

                quad_rule(w, t, x, pi, rule, prec);

                flag = f(y, t, param, order, prec);
                n_eval++;

                arb_mul(y, y, w, prec);
                arb_mul(y, y, h, prec);
                arb_add(w, result, y, prec);

                if (!arb_contains(w, result))
                {
                    arb_set(result, w);
                }
                else
                {
                    break;
                }

                if (pm > 0)
                {
                    arb_add(x, x, h, prec);
                }
                else
                {
                    arb_sub(x, x, h, prec);
                }

            }
        }

        arb_div_ui(x0, x0, 2, prec);

        if (j > 0)
        {
            arb_div_ui(h, h, 2, prec);
            arb_div_ui(result, result, 2, prec);
        }

        arb_sub(d + j, result, r + j, prec);
        arb_abs(d + j, d + j);

        if (arb_contains_zero(d + j))
            break;

        if (verbose > 0)
        {
            arb_log_base_ui(w, d + j, 10, prec);
            double d_1 = arf_get_d(arb_midref(w), ARF_RND_NEAR);
            if (j > 1)
            {
                arb_log_base_ui(w, d + j - 1, 10, prec);
                double d_0 = arf_get_d(arb_midref(w), ARF_RND_NEAR);
                flint_printf("\t%8.2f %8d %10.4f\n", d_1, n_eval, d_1 / d_0);
            }
            else
            {
                flint_printf("\t%8.2f %8d\n", d_1, n_eval);
            }
        }

        arb_set(r + j + 1, result);

    }

    arb_clear(x);
    arb_clear(y);
    arb_clear(h);
    arb_clear(x0);

    arb_clear(pi);

    arb_clear(w);
    arb_clear(t);

    _arb_vec_clear(r, max_iter);
    _arb_vec_clear(d, max_iter);
}

/*-- Example Integrands ------------------------------------------------------*/

int
f_rat(arb_ptr y, const arb_t x, void *param, slong order, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_pow_ui(y, x, 3, prec);
    arb_add_ui(y, y, 1, prec);
    arb_inv(y, y, prec);
    arb_clear(t);
    return 0;
}

void
f_rat_exc(arb_ptr exc, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_const_pi(exc, prec);
    arb_mul_ui(exc, exc, 2, prec);
    arb_set_ui(t, 27);
    arb_sqrt(t, t, prec);
    arb_div(exc, exc, t, prec);
    arb_clear(t);
}

int
g_rat(arb_ptr y, const arb_t x, void *param, slong order, slong prec)
{
    arb_pow_ui(y, x, 2, prec);
    arb_add_ui(y, y, 1, prec);
    arb_inv(y, y, prec);
    return 0;
}

void
g_rat_exc(arb_ptr exc, slong prec)
{
    arb_t t, w;
    arb_init(t);
    arb_init(w);
    arb_const_pi(exc, prec);
    arb_div_ui(exc, exc, 2, prec);
    arb_clear(t);
    arb_clear(w);
}

void
g_rat_exc_2(arb_ptr exc, slong prec)
{
    arb_const_pi(exc, prec);
}

int
f_log_exp(arb_ptr y, const arb_t x, void *param, slong order, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_neg(y, x);
    arb_exp(y, y, prec);
    arb_log(t, x, prec);
    arb_mul(y, y, t, prec);
    arb_clear(t);
    return 0;
}

void
f_log_exp_exc(arb_ptr exc, slong prec)
{
    arb_const_euler(exc, prec);
    arb_neg(exc, exc);
}

int
f_log_rat(arb_ptr y, const arb_t x, void *param, slong order, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_log(y, x, prec);
    arb_sqr(y, y, prec);
    arb_sqr(t, x, prec);
    arb_add_ui(t, t, 1, prec);
    arb_div(y, y, t, prec);
    arb_clear(t);
    return 0;
}

void
f_log_rat_exc(arb_ptr exc, slong prec)
{
    arb_const_pi(exc, prec);
    arb_pow_ui(exc, exc, 3, prec);
    arb_div_ui(exc, exc, 8, prec);
}

int
f_polylog(arb_ptr y, const arb_t x, void *param, slong order, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_sqr(y, x, prec);
    arb_exp(t, x, prec);
    arb_sub_ui(t, t, 1, prec);
    arb_div(y, y, t, prec);
    arb_clear(t);
    return 0;
}

void
f_polylog_exc(arb_ptr exc, slong prec)
{
    arb_zeta_ui(exc, 3, prec);
    arb_mul_ui(exc, exc, 2, prec);
}

int
f_exp_sqrt(arb_ptr y, const arb_t x, void *param, slong order, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_neg(y, x);
    arb_exp(y, y, prec);
    arb_sqrt(t, x, prec);
    arb_div(y, y, t, prec);
    arb_clear(t);
    return 0;
}

void
f_exp_sqrt_exc(arb_ptr exc, slong prec)
{
    arb_const_pi(exc, prec);
    arb_sqrt(exc, exc, prec);
}

int
f_gauss(arb_ptr y, const arb_t x, void *param, slong order, slong prec)
{
    arb_sqr(y, x, prec);
    arb_neg(y, y);
    arb_exp(y, y, prec);
    return 0;
}

void
f_gauss_exc(arb_ptr exc, slong prec)
{
    arb_const_pi(exc, prec);
    arb_sqrt(exc, exc, prec);
    arb_div_ui(exc, exc, 2, prec);
}

int
f_sqrt(arb_ptr y, const arb_t x, void *param, slong order, slong prec)
{
    arb_sqr(y, x, prec);
    arb_neg(y, y);
    arb_add_ui(y, y, 1, prec);
    arb_sqrt(y, y, prec);
    return 0;
}

void
f_sqrt_exc(arb_ptr exc, slong prec)
{
    arb_const_pi(exc, prec);
    arb_div_ui(exc, exc, 2, prec);
}

/*-- structure for integration examples --------------------------------------*/

typedef struct
{
    arb_calc_func_t fun;
    void (*exc)(arb_ptr out, slong prec);
    int rule;
    char *desc;
    char *result;
} example;

char *rule_desc[] = { "-1^+1", "0^oo", "-oo^+oo" };

example ex[] = {
    {f_rat, f_rat_exc, QUAD_SEMI_INF, "1/(1+x^3)", "2*pi/(3*sqrt(3))"},
    {f_log_exp, f_log_exp_exc, QUAD_SEMI_INF, "exp(-x)*log(x)", "-gamma"},
    {f_log_rat, f_log_rat_exc, QUAD_SEMI_INF, "log(x)^2/(1+x^2)", "pi^3/8"},
    {g_rat, g_rat_exc, QUAD_INTERVAL, "1/(1+x^2)", "pi/2"},
    {g_rat, g_rat_exc_2, QUAD_INF, "1/(1+x^2)", "pi"},
    {f_polylog, f_polylog_exc, QUAD_SEMI_INF, "x^2/(exp(x) - 1)", "2*zeta(3)"},
    {f_exp_sqrt, f_exp_sqrt_exc, QUAD_SEMI_INF, "exp(-x)/sqrt(x)", "sqrt(pi)"},
    {f_gauss, f_gauss_exc, QUAD_SEMI_INF, "exp(-x^2)", "sqrt(pi)/2"},
    {f_sqrt, f_sqrt_exc, QUAD_INTERVAL, "sqrt(1-x^2)", "pi/2"}
};

void
show_desc(int argc, char *argv[])
{
    flint_printf("Examples for double exponential integration.\n\n", argv[0]);
    flint_printf("usage: %s [-v] [-p precision] [-h]\n");
    flint_printf("\n\tcalculates a numerical approximation to:\n\n");
    for (slong j = 0; j < sizeof(ex) / sizeof(ex[0]); j++)
    {
        flint_printf("\tintegral_%s %s dx = %s\n",
                     rule_desc[ex[j].rule], ex[j].desc, ex[j].result);
    }
}

/*-- main program ------------------------------------------------------------*/

int
main(int argc, char *argv[])
{

    slong num_threads = 1;
    slong prec = 1024;

    int index, c, verbose = 0;

    opterr = 0;

    while ((c = getopt(argc, argv, "t:hvp:")) != -1)
    {
        switch (c)
        {
            case 'p':
                prec =
                    lround(strtol(optarg, (char **) NULL, 10) * log(10.) /
                           log(2.)) + 23;
                break;
            case 't':
                num_threads = atol(optarg);
                break;
            case 'v':
                verbose = 1;
                break;
            case 'h':
            default:
                show_desc(argc, argv);
                exit(0);
        }
    }

    if (num_threads > 1)
    {
        flint_set_num_threads(num_threads);
    }

    arb_t result, exc, d;

    arb_init(result);
    arb_init(exc);
    arb_init(d);

    for (slong j = 0; j < sizeof(ex) / sizeof(ex[0]); j++)
    {

        flint_printf("example: integral_%s %s dx = %s\n",
                     rule_desc[ex[j].rule], ex[j].desc, ex[j].result);

        quad(result, ex[j].fun, NULL, ex[j].rule, verbose, prec);

        ex[j].exc(exc, prec);

        arb_sub(d, result, exc, prec);
        arb_abs(d, d);

        flint_printf("\n\tresult = ");
        arb_printd(result, 16);
        flint_printf(", |result - exact| = ");
        arb_printn(d, 12, ARB_STR_NO_RADIUS);
        flint_printf("\n\n");

    }

}
