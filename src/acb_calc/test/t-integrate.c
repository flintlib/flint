/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_calc.h"

int
f_zero(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_zero(res);

    return 0;
}

int
f_indet(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_indeterminate(res);

    return 0;
}

int
f_cube(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_cube(res, z, prec);

    return 0;
}

int
f_sin(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_sin(res, z, prec);

    return 0;
}

int
f_sqrt(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_sqrt_analytic(res, z, order != 0, prec);

    return 0;
}

int
f_rsqrt(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_rsqrt_analytic(res, z, order != 0, prec);

    return 0;
}

int
f_log(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_log_analytic(res, z, order != 0, prec);

    return 0;
}

int
f_pow_pi(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_const_pi(res, prec);
    acb_pow_analytic(res, z, res, order != 0, prec);

    return 0;
}

int
f_circle(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_one(res);
    acb_submul(res, z, z, prec);
    acb_real_sqrtpos(res, res, order != 0, prec);

    return 0;
}

/* f(z) = sin(z + exp(z)) -- Rump's oscillatory example */
int
f_rump(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_exp(res, z, prec);
    acb_add(res, res, z, prec);
    acb_sin(res, res, prec);

    return 0;
}

/* f(z) = |z^4 + 10z^3 + 19z^2 - 6z - 6| exp(z)   (for real z) --
   Helfgott's integral on MathOverflow */
int
f_helfgott(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_add_si(res, z, 10, prec);
    acb_mul(res, res, z, prec);
    acb_add_si(res, res, 19, prec);
    acb_mul(res, res, z, prec);
    acb_add_si(res, res, -6, prec);
    acb_mul(res, res, z, prec);
    acb_add_si(res, res, -6, prec);

    acb_real_abs(res, res, order != 0, prec);

    if (acb_is_finite(res))
    {
        acb_t t;
        acb_init(t);
        acb_exp(t, z, prec);
        acb_mul(res, res, t, prec);
        acb_clear(t);
    }

    return 0;
}

/* f(z) = z sin(1/z), assume on real interval */
int
f_essing2(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    if ((order == 0) && acb_is_real(z) && arb_contains_zero(acb_realref(z)))
    {
        /* todo: arb_zero_pm_one, arb_unit_interval? */
        acb_zero(res);
        mag_one(arb_radref(acb_realref(res)));
    }
    else
    {
        acb_inv(res, z, prec);
        acb_sin(res, res, prec);
    }

    acb_mul(res, res, z, prec);

    return 0;
}

/* f(z) = sech(10(x-0.2))^2 + sech(100(x-0.4))^4 + sech(1000(x-0.6))^6 */
int
f_spike(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t a, b, c;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(a);
    acb_init(b);
    acb_init(c);

    acb_mul_ui(a, z, 10, prec);
    acb_sub_ui(a, a, 2, prec);
    acb_sech(a, a, prec);
    acb_pow_ui(a, a, 2, prec);

    acb_mul_ui(b, z, 100, prec);
    acb_sub_ui(b, b, 40, prec);
    acb_sech(b, b, prec);
    acb_pow_ui(b, b, 4, prec);

    acb_mul_ui(c, z, 1000, prec);
    acb_sub_ui(c, c, 600, prec);
    acb_sech(c, c, prec);
    acb_pow_ui(c, c, 6, prec);

    acb_add(res, a, b, prec);
    acb_add(res, res, c, prec);

    acb_clear(a);
    acb_clear(b);
    acb_clear(c);

    return 0;
}

/* more tests for the individual acb_real_* functions */

/* abs(sin(x))*cos(1+x) */
int
f_abs(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t s, c;
    acb_init(s);
    acb_init(c);
    acb_sin(s, z, prec);
    acb_add_ui(c, z, 1, prec);
    acb_cos(c, c, prec);
    acb_real_abs(s, s, order != 0, prec);
    acb_mul(res, s, c, prec);
    acb_clear(s);
    acb_clear(c);
    return 0;
}

/* sign(sin(x))*cos(1+x) */
int
f_sgn(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t s, c;
    acb_init(s);
    acb_init(c);
    acb_sin(s, z, prec);
    acb_add_ui(c, z, 1, prec);
    acb_cos(c, c, prec);
    acb_real_sgn(s, s, order != 0, prec);
    acb_mul(res, s, c, prec);
    acb_clear(s);
    acb_clear(c);
    return 0;
}

/* heaviside(sin(x))*cos(1+x) */
int
f_heaviside(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t s, c;
    acb_init(s);
    acb_init(c);
    acb_sin(s, z, prec);
    acb_add_ui(c, z, 1, prec);
    acb_cos(c, c, prec);
    acb_real_heaviside(s, s, order != 0, prec);
    acb_mul(res, s, c, prec);
    acb_clear(s);
    acb_clear(c);
    return 0;
}

/* floor(x-5) cos(1+x) */
int
f_floor(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t c;
    acb_init(c);
    acb_add_ui(c, z, 1, prec);
    acb_cos(c, c, prec);

    acb_sub_ui(res, z, 5, prec);
    acb_real_floor(res, res, order != 0, prec);

    acb_mul(res, res, c, prec);
    acb_clear(c);
    return 0;
}

/* ceil(x-5) cos(1+x) */
int
f_ceil(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t c;
    acb_init(c);
    acb_add_ui(c, z, 1, prec);
    acb_cos(c, c, prec);
    acb_sub_ui(res, z, 5, prec);
    acb_real_ceil(res, res, order != 0, prec);
    acb_mul(res, res, c, prec);
    acb_clear(c);
    return 0;
}

/* max(sin(x),cos(x)) */
int
f_max(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t s, c;
    acb_init(s);
    acb_init(c);
    acb_sin_cos(s, c, z, prec);
    acb_real_max(res, s, c, order != 0, prec);
    acb_clear(s);
    acb_clear(c);
    return 0;
}

/* min(sin(x),cos(x)) */
int
f_min(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t s, c;
    acb_init(s);
    acb_init(c);
    acb_sin_cos(s, c, z, prec);
    acb_real_min(res, s, c, order != 0, prec);
    acb_clear(s);
    acb_clear(c);
    return 0;
}

TEST_FUNCTION_START(acb_calc_integrate, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t ans, res, a, b, t;
        slong goal, prec, abs_goal;
        mag_t tol;
        acb_calc_integrate_opt_t opt;
        int integral;

        flint_set_num_threads(1 + n_randint(state, 3));

        acb_init(ans);
        acb_init(res);
        acb_init(a);
        acb_init(b);
        acb_init(t);
        mag_init(tol);
        acb_calc_integrate_opt_init(opt);

        goal = 2 + n_randint(state, 300);
        prec = 2 + n_randint(state, 300);
        abs_goal = n_randint(state, 300);
        mag_set_ui_2exp_si(tol, n_randint(state, 2), -abs_goal);

        opt->eval_limit = n_randint(state, 10000);

        if (n_randint(state, 2))
            opt->depth_limit = n_randint(state, 100);

        if (n_randint(state, 2))
            opt->deg_limit = n_randint(state, 100);

        opt->use_heap = n_randint(state, 2);

        integral = n_randint(state, 9);

        if (integral == 0)
        {
            acb_randtest(a, state, 1 + n_randint(state, 200), 2);
            acb_randtest(b, state, 1 + n_randint(state, 200), 2);
            acb_cos(ans, a, prec);
            acb_cos(res, b, prec);
            acb_sub(ans, ans, res, prec);
            acb_calc_integrate(res, f_sin, NULL, a, b, goal, tol, opt, prec);
        }
        else if (integral == 1)
        {
            acb_randtest(a, state, 1 + n_randint(state, 200), 2);
            acb_randtest(b, state, 1 + n_randint(state, 200), 2);
            acb_pow_ui(ans, a, 4, prec);
            acb_pow_ui(res, b, 4, prec);
            acb_sub(ans, res, ans, prec);
            acb_mul_2exp_si(ans, ans, -2);
            acb_calc_integrate(res, f_cube, NULL, a, b, goal, tol, opt, prec);
        }
        else if (integral == 2)
        {
            acb_zero(a);
            acb_one(b);
            acb_const_pi(ans, prec);
            acb_mul_2exp_si(ans, ans, -2);
            if (n_randint(state, 2))
            {
                acb_swap(a, b);
                acb_neg(ans, ans);
            }
            acb_calc_integrate(res, f_circle, NULL, a, b, goal, tol, opt, prec);
        }
        else if (integral == 3)
        {
            acb_zero(a);
            acb_one(b);
            acb_zero(ans);
            arb_set_str(acb_realref(ans), "0.3785300171241613098817352756283519095343133 +/- 1e-40", prec);
            if (n_randint(state, 2))
            {
                acb_swap(a, b);
                acb_neg(ans, ans);
            }
            acb_calc_integrate(res, f_essing2, NULL, a, b, goal, tol, opt, prec);
        }
        else if (integral == 4)
        {
            acb_zero(a);
            acb_one(b);
            acb_zero(ans);
            arb_set_str(acb_realref(ans), "11.147310550057139733915902084255301415775813549800589418261584268232061665808482234384871404010464 +/- 2.28e-97", prec);
            if (n_randint(state, 2))
            {
                acb_swap(a, b);
                acb_neg(ans, ans);
            }
            acb_calc_integrate(res, f_helfgott, NULL, a, b, goal, tol, opt, prec);
        }
        else if (integral == 5)
        {
            acb_zero(a);
            acb_set_ui(b, 8);
            acb_zero(ans);
            arb_set_str(acb_realref(ans), "0.34740017265724780787951215911989312465745625486618018388549271361674821398878532052968510434660 +/- 5.97e-96", prec);
            if (n_randint(state, 2))
            {
                acb_swap(a, b);
                acb_neg(ans, ans);
            }
            acb_calc_integrate(res, f_rump, NULL, a, b, goal, tol, opt, prec);
        }
        else if (integral == 6)
        {
            acb_zero(a);
            acb_one(b);
            acb_zero(ans);
            arb_set_str(acb_realref(ans), "0.21080273550054927737564325570572915436090918643678119034785050587872061312814550020505868926155764 +/- 3.72e-99", prec);
            if (n_randint(state, 2))
            {
                acb_swap(a, b);
                acb_neg(ans, ans);
            }
            acb_calc_integrate(res, f_spike, NULL, a, b, goal, tol, opt, prec);
        }
        else if (integral == 7)
        {
            acb_randtest_special(a, state, 1 + n_randint(state, 200), 2);
            acb_randtest_special(b, state, 1 + n_randint(state, 200), 2);
            acb_zero(ans);
            acb_calc_integrate(res, f_zero, NULL, a, b, goal, tol, opt, prec);
        }
        else if (integral == 8)
        {
            acb_randtest_special(a, state, 1 + n_randint(state, 200), 2);
            acb_randtest_special(b, state, 1 + n_randint(state, 200), 2);
            acb_indeterminate(ans);
            acb_calc_integrate(res, f_indet, NULL, a, b, goal, tol, opt, prec);
        }

        if (!acb_overlaps(res, ans))
        {
            flint_printf("FAIL! (iter = %wd)\n", iter);
            flint_printf("integral = %d, prec = %wd, goal = %wd\n", integral, prec, goal);
            flint_printf("a = "); acb_printn(a, 150, 0); flint_printf("\n");
            flint_printf("b = "); acb_printn(b, 150, 0); flint_printf("\n");
            flint_printf("res = "); acb_printn(res, 150, 0); flint_printf("\n\n");
            flint_printf("ans = "); acb_printn(ans, 150, 0); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(ans);
        acb_clear(res);
        acb_clear(a);
        acb_clear(b);
        acb_clear(t);
        mag_clear(tol);
    }

    /* more tests for the individual real extensions and branched functions */
    {
        acb_t a, b, z, w;
        slong prec;
        mag_t tol;

        acb_init(a);
        acb_init(b);
        acb_init(z);
        acb_init(w);
        mag_init(tol);

        acb_zero(a);
        acb_set_ui(b, 10);

        prec = 53;
        mag_set_ui_2exp_si(tol, 1, -prec);

        acb_calc_integrate(z, f_abs, NULL, a, b, prec, tol, NULL, prec);
        arb_set_str(acb_realref(w), "-1.3517710956465318592 +/- 1e-17", prec);

        if (!acb_overlaps(z, w) || acb_rel_accuracy_bits(z) < prec - 10)
        {
            flint_printf("FAIL (abs)\n");
            flint_printf("z = "); acb_printn(z, 20,  0); flint_printf("\n");
            flint_printf("w = "); acb_printn(w, 20,  0); flint_printf("\n");
            flint_abort();
        }

        acb_calc_integrate(z, f_sgn, NULL, a, b, prec, tol, NULL, prec);
        arb_set_str(acb_realref(w), "-4.8903066871045720895 +/- 1e-17", prec);

        if (!acb_overlaps(z, w) || acb_rel_accuracy_bits(z) < prec - 10)
        {
            flint_printf("FAIL (sgn)\n");
            flint_printf("z = "); acb_printn(z, 20,  0); flint_printf("\n");
            flint_printf("w = "); acb_printn(w, 20,  0); flint_printf("\n");
            flint_abort();
        }

        acb_calc_integrate(z, f_heaviside, NULL, a, b, prec, tol, NULL, prec);
        arb_set_str(acb_realref(w), "-3.3658839392315860266 +/- 1e-17", prec);

        if (!acb_overlaps(z, w) || acb_rel_accuracy_bits(z) < prec - 10)
        {
            flint_printf("FAIL (heaviside)\n");
            flint_printf("z = "); acb_printn(z, 20,  0); flint_printf("\n");
            flint_printf("w = "); acb_printn(w, 20,  0); flint_printf("\n");
            flint_abort();
        }

        acb_calc_integrate(z, f_floor, NULL, a, b, prec, tol, NULL, prec);
        arb_set_str(acb_realref(w), "-0.36232328857344524392 +/- 1e-17", prec);

        if (!acb_overlaps(z, w) || acb_rel_accuracy_bits(z) < prec - 16)
        {
            flint_printf("FAIL (floor)\n");
            flint_printf("z = "); acb_printn(z, 20,  0); flint_printf("\n");
            flint_printf("w = "); acb_printn(w, 20,  0); flint_printf("\n");
            flint_abort();
        }

        acb_calc_integrate(z, f_ceil, NULL, a, b, prec, tol, NULL, prec);
        arb_set_str(acb_realref(w), "-2.2037844799320452076 +/- 1e-17", prec);

        if (!acb_overlaps(z, w) || acb_rel_accuracy_bits(z) < prec - 16)
        {
            flint_printf("FAIL (ceil)\n");
            flint_printf("z = "); acb_printn(z, 20,  0); flint_printf("\n");
            flint_printf("w = "); acb_printn(w, 20,  0); flint_printf("\n");
            flint_abort();
        }

        acb_calc_integrate(z, f_max, NULL, a, b, prec, tol, NULL, prec);
        arb_set_str(acb_realref(w), "5.0817122161957375987 +/- 1e-17", prec);

        if (!acb_overlaps(z, w) || acb_rel_accuracy_bits(z) < prec - 10)
        {
            flint_printf("FAIL (max)\n");
            flint_printf("z = "); acb_printn(z, 20,  0); flint_printf("\n");
            flint_printf("w = "); acb_printn(w, 20,  0); flint_printf("\n");
            flint_abort();
        }

        acb_calc_integrate(z, f_min, NULL, a, b, prec, tol, NULL, prec);
        arb_set_str(acb_realref(w), "-3.7866617980086549598 +/- 1e-17", prec);

        if (!acb_overlaps(z, w) || acb_rel_accuracy_bits(z) < prec - 10)
        {
            flint_printf("FAIL (min)\n");
            flint_printf("z = "); acb_printn(z, 20,  0); flint_printf("\n");
            flint_printf("w = "); acb_printn(w, 20,  0); flint_printf("\n");
            flint_abort();
        }

        prec = 64;
        mag_set_ui_2exp_si(tol, 1, -prec);

        acb_set_ui(a, 10);
        acb_set_ui(b, 11);
        acb_calc_integrate(z, f_sqrt, NULL, a, b, prec, tol, NULL, prec);
        arb_set_str(acb_realref(w), "3.2400640614837366802 +/- 1.65e-20", prec);

        if (!acb_overlaps(z, w) || acb_rel_accuracy_bits(z) < prec - 10)
        {
            flint_printf("FAIL (sqrt)\n");
            flint_printf("z = "); acb_printn(z, 20,  0); flint_printf("\n");
            flint_printf("w = "); acb_printn(w, 20,  0); flint_printf("\n");
            flint_abort();
        }

        acb_calc_integrate(z, f_rsqrt, NULL, a, b, prec, tol, NULL, prec);
        arb_set_str(acb_realref(w), "0.30869426037404103423 +/- 2.08e-21", prec);

        if (!acb_overlaps(z, w) || acb_rel_accuracy_bits(z) < prec - 10)
        {
            flint_printf("FAIL (rsqrt)\n");
            flint_printf("z = "); acb_printn(z, 20,  0); flint_printf("\n");
            flint_printf("w = "); acb_printn(w, 20,  0); flint_printf("\n");
            flint_abort();
        }

        acb_calc_integrate(z, f_log, NULL, a, b, prec, tol, NULL, prec);
        arb_set_str(acb_realref(w), "2.3509970708416191445 +/- 1.47e-21", prec);

        if (!acb_overlaps(z, w) || acb_rel_accuracy_bits(z) < prec - 10)
        {
            flint_printf("FAIL (log)\n");
            flint_printf("z = "); acb_printn(z, 20,  0); flint_printf("\n");
            flint_printf("w = "); acb_printn(w, 20,  0); flint_printf("\n");
            flint_abort();
        }

        acb_calc_integrate(z, f_pow_pi, NULL, a, b, prec, tol, NULL, prec);
        arb_set_str(acb_realref(w), "1619.0628341842463436 +/- 2.79e-17", prec);

        if (!acb_overlaps(z, w) || acb_rel_accuracy_bits(z) < prec - 10)
        {
            flint_printf("FAIL (pow_pi)\n");
            flint_printf("z = "); acb_printn(z, 20,  0); flint_printf("\n");
            flint_printf("w = "); acb_printn(w, 20,  0); flint_printf("\n");
            flint_abort();
        }

        acb_set_d_d(a, -10, -1);
        acb_set_d_d(b, -10, +1);
        acb_zero(w);

        acb_calc_integrate(z, f_sqrt, NULL, a, b, prec, tol, NULL, prec);
        arb_set_str(acb_imagref(w), "0.15801534879296271761 +/- 4.91e-21", prec);

        if (!acb_overlaps(z, w) || acb_rel_accuracy_bits(z) < prec - 25)
        {
            flint_printf("FAIL (sqrt 2)\n");
            flint_printf("z = "); acb_printn(z, 20,  0); flint_printf("\n");
            flint_printf("w = "); acb_printn(w, 20,  0); flint_printf("\n");
            flint_abort();
        }

        acb_calc_integrate(z, f_rsqrt, NULL, a, b, prec, tol, NULL, prec);
        arb_set_str(acb_imagref(w), "0.015762235473606564294 +/- 2.80e-22", prec);

        if (!acb_overlaps(z, w) || acb_rel_accuracy_bits(z) < prec - 25)
        {
            flint_printf("FAIL (rsqrt 2)\n");
            flint_printf("z = "); acb_printn(z, 20,  0); flint_printf("\n");
            flint_printf("w = "); acb_printn(w, 20,  0); flint_printf("\n");
            flint_abort();
        }

        acb_calc_integrate(z, f_log, NULL, a, b, prec, tol, NULL, prec);
        arb_set_str(acb_imagref(w), "4.6084935666644999985 +/- 4.69e-20", prec);

        if (!acb_overlaps(z, w) || acb_rel_accuracy_bits(z) < prec - 25)
        {
            flint_printf("FAIL (log 2)\n");
            flint_printf("z = "); acb_printn(z, 20,  0); flint_printf("\n");
            flint_printf("w = "); acb_printn(w, 20,  0); flint_printf("\n");
            flint_abort();
        }

        acb_calc_integrate(z, f_pow_pi, NULL, a, b, prec, tol, NULL, prec);
        arb_set_str(acb_imagref(w), "-2660.1245860252382407 +/- 8.64e-18", prec);

        if (!acb_overlaps(z, w) || acb_rel_accuracy_bits(z) < prec - 25)
        {
            flint_printf("FAIL (pow_pi 2)\n");
            flint_printf("z = "); acb_printn(z, 20,  0); flint_printf("\n");
            flint_printf("w = "); acb_printn(w, 20,  0); flint_printf("\n");
            flint_abort();
        }

        prec = 100;
        mag_set_ui_2exp_si(tol, 1, -prec);

        acb_set_d_d(a, -100, -1);
        acb_set_d_d(b, -100, +2);
        acb_zero(w);
        acb_calc_integrate(z, f_log, NULL, a, b, prec, tol, NULL, prec);
        arb_set_str(acb_realref(w), "-3.12659390337983876281336380373 +/- 3.97e-30", prec);
        arb_set_str(acb_imagref(w), "13.8156605414673448203655975232 +/- 2.39e-29", prec);

        if (!acb_overlaps(z, w) || acb_rel_accuracy_bits(z) < prec - 25)
        {
            flint_printf("FAIL (log 3)\n");
            flint_printf("z = "); acb_printn(z, 20,  0); flint_printf("\n");
            flint_printf("w = "); acb_printn(w, 20,  0); flint_printf("\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(z);
        acb_clear(w);
        mag_clear(tol);
    }

    TEST_FUNCTION_END(state);
}
