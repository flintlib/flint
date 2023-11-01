/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <float.h>
#include <mpfr.h>
#include "ulong_extras.h"
#include "double_extras.h"

#define ONE_OVER_E ldexp(6627126856707895.0, -54)

TEST_FUNCTION_START(d_lambertw, state)
{
    double x, w, tol;
    slong iter, prec = 70;
    mpfr_t xx, ww, wnew, t, u, v, p, q, max_err;

    mpfr_init2(xx, prec);
    mpfr_init2(ww, prec);
    mpfr_init2(wnew, prec);
    mpfr_init2(t, prec);
    mpfr_init2(u, prec);
    mpfr_init2(v, prec);
    mpfr_init2(p, prec);
    mpfr_init2(q, prec);
    mpfr_init2(max_err, prec);
    mpfr_set_ui(max_err, 0, MPFR_RNDN);

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        x = d_randtest(state);

        switch (n_randint(state, 3))
        {
            /* singularity near -1/e */
            case 0:
                x = ldexp(x, -n_randint(state, -DBL_MIN_EXP+1));
                x = -ONE_OVER_E + x;
                tol = 50 * DBL_EPSILON;
                break;
            /* negative, not close to -1/e */
            case 1:
                x = d_randtest(state);
                x = ldexp(x, -n_randint(state, -DBL_MIN_EXP+1));
                x = x * -(1./4);
                tol = 2 * DBL_EPSILON;
                break;
            /* positive */
            default:
                x = d_randtest(state);
                x = ldexp(x, (int) n_randint(state, DBL_MAX_EXP-DBL_MIN_EXP-1)
                        + DBL_MIN_EXP);
                tol = 2 * DBL_EPSILON;
                break;
        }

        w = d_lambertw(x);

        mpfr_set_d(xx, x, MPFR_RNDN);
        mpfr_set_d(ww, w, MPFR_RNDN);

        /* t = exp(w) */
        mpfr_exp(t, ww, MPFR_RNDN);

        /* u = 2*w + 2 */
        mpfr_mul_ui(u, ww, 2, MPFR_RNDN);
        mpfr_add_ui(u, u, 2, MPFR_RNDN);

        /* v = w*t - x */
        mpfr_mul(v, t, ww, MPFR_RNDN);
        mpfr_sub(v, v, xx, MPFR_RNDN);

        /* p = u * v */
        mpfr_mul(p, u, v, MPFR_RNDN);

        /* q = (u*t*(w+1) - (w+2)*v) */
        mpfr_mul(q, u, t, MPFR_RNDN);
        mpfr_add_ui(t, ww, 1, MPFR_RNDN);
        mpfr_mul(q, q, t, MPFR_RNDN);

        mpfr_add_ui(t, ww, 2, MPFR_RNDN);
        mpfr_mul(t, t, v, MPFR_RNDN);
        mpfr_sub(q, q, t, MPFR_RNDN);

        /* wnew = w - p / q */
        mpfr_div(p, p, q, MPFR_RNDN);
        mpfr_sub(wnew, ww, p, MPFR_RNDN);

        /* relative error */
        mpfr_sub(t, ww, wnew, MPFR_RNDA);
        mpfr_div(t, t, wnew, MPFR_RNDA);
        mpfr_abs(t, t, MPFR_RNDA);

        if (mpfr_get_d(t, MPFR_RNDA) > tol)
        {
            flint_printf("FAIL\n");
            flint_printf("x = %.17g, w = %.17g, error = %g\n", x, w,
                mpfr_get_d(t, MPFR_RNDA));
            fflush(stdout);
            flint_abort();
        }

#if 0
        if (mpfr_cmp(t, max_err) > 0)
        {
            flint_printf("new record: ");
            flint_printf("x=%.20g w=%.20g wnew=%.20g relative error: %g\n",
                x, w, mpfr_get_d(wnew, MPFR_RNDN), mpfr_get_d(t, MPFR_RNDN));
            mpfr_set(max_err, t, MPFR_RNDN);
        }
#endif

    }

    mpfr_clear(xx);
    mpfr_clear(ww);
    mpfr_clear(wnew);
    mpfr_clear(t);
    mpfr_clear(u);
    mpfr_clear(v);
    mpfr_clear(p);
    mpfr_clear(q);
    mpfr_clear(max_err);

    mpfr_free_cache();

    TEST_FUNCTION_END(state);
}
