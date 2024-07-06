/*
    Copyright (C) 2024 Matthias Gessinger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

void _fmpz_sum_powers_naive_horner(fmpz_t f, const fmpz_t g, ulong exp)
{
	fmpz_t t;
	fmpz_init(t);

	fmpz_zero(f);

	for (ulong i = 0; i <= exp; i++)
	{
		fmpz_pow_ui(t, g, i);
		fmpz_add(f, f, t);
	}

	fmpz_clear(t);
}

TEST_FUNCTION_START(fmpz_sum_powers_horner, state)
{
	int i, result;

	for (i = 0; i < 1000 * flint_test_multiplier(); i++)
	{
		fmpz_t a, b, c;
		ulong exp;
		int aliasing;

		fmpz_init(a);
		fmpz_init(b);
		fmpz_init(c);

		fmpz_randtest(a, state, 100);

		exp = n_randint(state, 250);

		aliasing = n_randint(state, 2);

		/* The reference result */
		_fmpz_sum_powers_naive_horner(c, a, exp);

		/* Horner-style evaluation */
		if (aliasing == 1)
		{
			fmpz_sum_powers_horner(a, a, exp);
			fmpz_set(b, a);
		}
		else
		{
			fmpz_sum_powers_horner(b, a, exp);
		}

		result = (fmpz_cmp(c, b) == 0);

		if (!result)
		{
			flint_printf("FAIL in horner evaluation:\n");

			flint_printf("Expected: ");
			fmpz_print(c);
			flint_printf("\n");

			flint_printf("Actual: ");
			fmpz_print(b);
			flint_printf("\n");

			fflush(stdout);
			flint_abort();
		}

		fmpz_clear(a);
		fmpz_clear(b);
		fmpz_clear(c);
	}

	TEST_FUNCTION_END(state);
}