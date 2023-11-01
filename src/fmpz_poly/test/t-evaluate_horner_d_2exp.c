/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpfr.h"
#include "double_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_evaluate_horner_d_2exp, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f;
        double x, y, z, t;
        slong xexp, yexp;

        x = d_randtest(state);
        xexp = n_randint(state, 20) - 10;

        fmpz_poly_init(f);
        fmpz_poly_randtest(f, state, 1 + n_randint(state, 40), 1 + n_randint(state, 100));
        fmpz_poly_scalar_abs(f, f);

        y = fmpz_poly_evaluate_horner_d_2exp2(&yexp, f, x, xexp);
        z = fmpz_poly_evaluate_horner_d(f, ldexp(x, xexp));
        t = ldexp(y, yexp);

        if (fabs(t - z) > 1e-13 * fabs(z))
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(f), flint_printf("\n\n");
            flint_printf("x, xexp = %.20g  %wd\n\n", x, xexp);
            flint_printf("y, yexp = %.20g  %wd\n\n", y, yexp);
            flint_printf("z = %.20g\n\n", z);
            flint_printf("y = %.20g\n\n", t);
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f;
        double x, y;
        slong xexp, yexp, i;
        mpfr_t z, s, t, u, v, w, e;

        x = d_randtest(state);
        xexp = n_randint(state, 2000) - 1000;

        fmpz_poly_init(f);
        fmpz_poly_randtest(f, state, 1 + n_randint(state, 100), 1 + n_randint(state, 1000));
        fmpz_poly_scalar_abs(f, f);
        mpfr_init2(z, 64);
        mpfr_init2(s, 64);
        mpfr_init2(t, 64);
        mpfr_init2(u, 64);
        mpfr_init2(v, 64);
        mpfr_init2(w, 64);
        mpfr_init2(e, 64);

        y = fmpz_poly_evaluate_horner_d_2exp2(&yexp, f, x, xexp);

        mpfr_set_d(z, x, MPFR_RNDN);
        mpfr_mul_2si(z, z, xexp, MPFR_RNDN);
        mpfr_set_ui(s, 0, MPFR_RNDN);
        mpfr_set_ui(t, 1, MPFR_RNDN);

        for (i = 0; i < f->length; i++)
        {
            fmpz_get_mpfr(u, f->coeffs + i, MPFR_RNDN);
            mpfr_mul(u, u, t, MPFR_RNDN);
            mpfr_add(s, s, u, MPFR_RNDN);
            mpfr_mul(t, t, z, MPFR_RNDN);
        }

        mpfr_set_d(v, y, MPFR_RNDN);
        mpfr_mul_2si(v, v, yexp, MPFR_RNDN);

        mpfr_sub(e, s, v, MPFR_RNDN);
        mpfr_abs(e, e, MPFR_RNDN);

        mpfr_abs(w, s, MPFR_RNDN);
        mpfr_mul_ui(w, w, f->length + 1, MPFR_RNDN);
        mpfr_mul_2si(w, w, -51, MPFR_RNDN);

        if (mpfr_cmp(e, w) > 0)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(f), flint_printf("\n\n");
            mpfr_printf("%.17Rg\n", s);
            mpfr_printf("%.17Rg\n", v);
            mpfr_printf("%.17Rg\n", e);
            mpfr_printf("%.17Rg\n\n", w);
            fflush(stdout);
            flint_abort();
        }

        mpfr_clear(z);
        mpfr_clear(s);
        mpfr_clear(t);
        mpfr_clear(u);
        mpfr_clear(v);
        mpfr_clear(w);
        mpfr_clear(e);
        fmpz_poly_clear(f);
    }

    TEST_FUNCTION_END(state);
}
