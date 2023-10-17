/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "fmpz.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_exp_series, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* Check aliasing of a and c */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        slong n = n_randint(state, 50) + 1;

        fmpq_poly_init(a);
        fmpq_poly_init(b);

        fmpq_poly_randtest_not_zero(a, state, n_randint(state, 50) + 1, 50);
        fmpq_poly_set_coeff_ui(a, 0, UWORD(0));

        fmpq_poly_canonicalise(a);

        fmpq_poly_exp_series(b, a, n);
        fmpq_poly_exp_series(a, a, n);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(a, b) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Check exp(a+b) = exp(a) * exp(b) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, ab, expa, expb, expab, expa_expb;
        slong n = n_randint(state, 20) + 1;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(ab);
        fmpq_poly_init(expa);
        fmpq_poly_init(expb);
        fmpq_poly_init(expab);
        fmpq_poly_init(expa_expb);

        fmpq_poly_randtest(a, state, n_randint(state, 20) + 1, 1 + n_randint(state, 80));
        fmpq_poly_set_coeff_ui(a, 0, UWORD(0));

        fmpq_poly_randtest(b, state, n_randint(state, 20) + 1, 1 + n_randint(state, 80));
        fmpq_poly_set_coeff_ui(b, 0, UWORD(0));

        fmpq_poly_add(ab, a, b);

        fmpq_poly_exp_series(expab, ab, n);
        fmpq_poly_exp_series(expa, a, n);
        fmpq_poly_exp_series(expb, b, n);
        fmpq_poly_mullow(expa_expb, expa, expb, n);

        cflags |= fmpq_poly_is_canonical(expa) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(expb) ? 0 : 2;
        cflags |= fmpq_poly_is_canonical(expab) ? 0 : 4;
        result = (fmpq_poly_equal(expab, expa_expb) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n\n");
            flint_printf("b = "), fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("exp(a) = "), fmpq_poly_debug(expa), flint_printf("\n\n");
            flint_printf("exp(b) = "), fmpq_poly_debug(expb), flint_printf("\n\n");
            flint_printf("exp(ab) = "), fmpq_poly_debug(expab), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(ab);
        fmpq_poly_clear(expa);
        fmpq_poly_clear(expb);
        fmpq_poly_clear(expab);
        fmpq_poly_clear(expa_expb);
    }

    /* Compare with nmod_poly_exp */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, expa;
        nmod_poly_t pa, pexpa, pexpa2;
        ulong p;
        slong n = n_randint(state, 200) + 1;

        fmpq_poly_init(a);
        fmpq_poly_init(expa);

        p = n_nextprime((UWORD(1) << (FLINT_BITS - 1)) + 1000, 0);

        nmod_poly_init(pa, p);
        nmod_poly_init(pexpa, p);
        nmod_poly_init(pexpa2, p);

        fmpq_poly_randtest(a, state, n_randint(state, 200) + 1, 1 + n_randint(state, 200));
        fmpz_randtest_not_zero(fmpq_poly_denref(a), state, 1 + n_randint(state, 200));
        fmpq_poly_canonicalise(a);

        fmpq_poly_set_coeff_ui(a, 0, UWORD(0));

        fmpq_poly_get_nmod_poly(pa, a);
        nmod_poly_exp_series(pexpa, pa, n);

        if (n_randint(state, 2))
        {
            fmpq_poly_exp_series(expa, a, n);
        }
        else
        {
            fmpq_poly_set(expa, a);
            fmpq_poly_exp_series(expa, expa, n);
        }

        fmpq_poly_get_nmod_poly(pexpa2, expa);

        if (!fmpq_poly_is_canonical(expa) || !nmod_poly_equal(pexpa, pexpa2))
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n\n");
            flint_printf("exp(a) = "), fmpq_poly_debug(expa), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(expa);
        nmod_poly_clear(pa);
        nmod_poly_clear(pexpa);
        nmod_poly_clear(pexpa2);
    }

    TEST_FUNCTION_END(state);
}
