/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "arb_calc.h"
#include "acb_calc.h"

#ifdef __GNUC__
# define log __builtin_log
#else
# include <math.h>
#endif

int
acb_calc_integrate_taylor(acb_t res,
    acb_calc_func_t func, void * param,
    const acb_t a, const acb_t b,
    const arf_t inner_radius,
    const arf_t outer_radius,
    slong accuracy_goal, slong prec)
{
    slong num_steps, step, N, bp;
    int result;

    acb_t delta, m, x, y1, y2, sum;
    acb_ptr taylor_poly;
    arf_t err;

    acb_init(delta);
    acb_init(m);
    acb_init(x);
    acb_init(y1);
    acb_init(y2);
    acb_init(sum);
    arf_init(err);

    acb_sub(delta, b, a, prec);

    /* precision used for bounds calculations */
    bp = MAG_BITS;

    /* compute the number of steps */
    {
        arf_t t;
        arf_init(t);
        acb_get_abs_ubound_arf(t, delta, bp);
        arf_div(t, t, inner_radius, bp, ARF_RND_UP);
        arf_mul_2exp_si(t, t, -1);
        num_steps = (slong) (arf_get_d(t, ARF_RND_UP) + 1.0);
        /* make sure it's not something absurd */
        num_steps = FLINT_MIN(num_steps, 10 * prec);
        num_steps = FLINT_MAX(num_steps, 1);
        arf_clear(t);
    }

    result = ARB_CALC_SUCCESS;

    acb_zero(sum);

    for (step = 0; step < num_steps; step++)
    {
        /* midpoint of subinterval */
        acb_mul_ui(m, delta, 2 * step + 1, prec);
        acb_div_ui(m, m, 2 * num_steps, prec);
        acb_add(m, m, a, prec);

        if (arb_calc_verbose)
        {
            flint_printf("integration point %wd/%wd: ", 2 * step + 1, 2 * num_steps);
            acb_printd(m, 15); flint_printf("\n");
        }

        /* evaluate at +/- x */
        /* TODO: exactify m, and include error in x? */
        acb_div_ui(x, delta, 2 * num_steps, prec);

        /* compute bounds and number of terms to use */
        {
            arb_t cbound, xbound, rbound;
            arf_t C, D, R, X, T;
            double DD, TT, NN;

            arb_init(cbound);
            arb_init(xbound);
            arb_init(rbound);
            arf_init(C);
            arf_init(D);
            arf_init(R);
            arf_init(X);
            arf_init(T);

            /* R is the outer radius */
            arf_set(R, outer_radius);

            /* X = upper bound for |x| */
            acb_get_abs_ubound_arf(X, x, bp);
            arb_set_arf(xbound, X);

            /* Compute C(m,R). Important subtlety: due to rounding when
               computing m, we will in general be farther than R away from
               the integration path. But since acb_calc_cauchy_bound
               actually integrates over the area traced by a complex
               interval, it will catch any extra singularities (giving
               an infinite bound). */
            arb_set_arf(rbound, outer_radius);
            acb_calc_cauchy_bound(cbound, func, param, m, rbound, 8, bp);
            arf_set_mag(C, arb_radref(cbound));
            arf_add(C, arb_midref(cbound), C, bp, ARF_RND_UP);

            /* Sanity check: we need C < inf and R > X */
            if (arf_is_finite(C) && arf_cmp(R, X) > 0)
            {
                /* Compute upper bound for D = C * R * X / (R - X) */
                arf_mul(D, C, R, bp, ARF_RND_UP);
                arf_mul(D, D, X, bp, ARF_RND_UP);
                arf_sub(T, R, X, bp, ARF_RND_DOWN);
                arf_div(D, D, T, bp, ARF_RND_UP);

                /* Compute upper bound for T = (X / R) */
                arf_div(T, X, R, bp, ARF_RND_UP);

                /* Choose N */
                /* TODO: use arf arithmetic to avoid overflow */
                /* TODO: use relative accuracy (look at |f(m)|?) */
                DD = arf_get_d(D, ARF_RND_UP);
                TT = arf_get_d(T, ARF_RND_UP);
                NN = -(accuracy_goal * 0.69314718055994530942 + log(DD)) / log(TT);
                N = NN + 0.5;
                N = FLINT_MIN(N, 100 * prec);
                N = FLINT_MAX(N, 1);

                /* Tail bound: D / (N + 1) * T^N */
                {
                    mag_t TT;
                    mag_init(TT);
                    arf_get_mag(TT, T);
                    mag_pow_ui(TT, TT, N);
                    arf_set_mag(T, TT);
                    mag_clear(TT);
                }
                arf_mul(D, D, T, bp, ARF_RND_UP);
                arf_div_ui(err, D, N + 1, bp, ARF_RND_UP);
            }
            else
            {
                N = 1;
                arf_pos_inf(err);
                result = ARB_CALC_NO_CONVERGENCE;
            }

            if (arb_calc_verbose)
            {
                flint_printf("N = %wd; bound: ", N); arf_printd(err, 15); flint_printf("\n");
                flint_printf("R: "); arf_printd(R, 15); flint_printf("\n");
                flint_printf("C: "); arf_printd(C, 15); flint_printf("\n");
                flint_printf("X: "); arf_printd(X, 15); flint_printf("\n");
            }

            arb_clear(cbound);
            arb_clear(xbound);
            arb_clear(rbound);
            arf_clear(C);
            arf_clear(D);
            arf_clear(R);
            arf_clear(X);
            arf_clear(T);
        }

        /* evaluate Taylor polynomial */
        taylor_poly = _acb_vec_init(N + 1);
        func(taylor_poly, m, param, N, prec);
        _acb_poly_integral(taylor_poly, taylor_poly, N + 1, prec);
        _acb_poly_evaluate(y2, taylor_poly, N + 1, x, prec);
        acb_neg(x, x);
        _acb_poly_evaluate(y1, taylor_poly, N + 1, x, prec);
        acb_neg(x, x);

        /* add truncation error */
        arb_add_error_arf(acb_realref(y1), err);
        arb_add_error_arf(acb_imagref(y1), err);
        arb_add_error_arf(acb_realref(y2), err);
        arb_add_error_arf(acb_imagref(y2), err);

        acb_add(sum, sum, y2, prec);
        acb_sub(sum, sum, y1, prec);

        if (arb_calc_verbose)
        {
            flint_printf("values:  ");
            acb_printd(y1, 15); flint_printf("  ");
            acb_printd(y2, 15); flint_printf("\n");
        }

        _acb_vec_clear(taylor_poly, N + 1);

        if (result == ARB_CALC_NO_CONVERGENCE)
            break;
    }

    acb_set(res, sum);

    acb_clear(delta);
    acb_clear(m);
    acb_clear(x);
    acb_clear(y1);
    acb_clear(y2);
    acb_clear(sum);
    arf_clear(err);

    return result;
}

