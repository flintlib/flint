/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_xgcd_modular, state)
{
    int i, result;

    /* Check s * f + t * g == r */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t d, f, g, s, t;
        fmpz_poly_t tmp, res;
        fmpz_t r;
        int type;

        fmpz_poly_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(s);
        fmpz_poly_init(t);
        fmpz_poly_init(tmp);
        fmpz_poly_init(res);
        fmpz_init(r);

        do
        {
            fmpz_poly_randtest(f, state, n_randint(state, 50), 100);
            fmpz_poly_randtest(g, state, n_randint(state, 50), 100);
            fmpz_poly_primitive_part(f, f);
            fmpz_poly_primitive_part(g, g);
            fmpz_poly_gcd(d, f, g);
        } while (fmpz_poly_length(d) != 1);

        type = n_randint(state, 5);

        switch (type)
        {
            case 0: /* No aliasing */
                fmpz_poly_xgcd_modular(r, s, t, f, g);
                break;

            case 1: /* s aliased with f */
                fmpz_poly_set(s, f);
                fmpz_poly_xgcd_modular(r, s, t, s, g);
                break;

            case 2: /* s aliased with g */
                fmpz_poly_set(s, g);
                fmpz_poly_xgcd_modular(r, s, t, f, s);
                break;

            case 3: /* t aliased with f */
                fmpz_poly_set(t, f);
                fmpz_poly_xgcd_modular(r, s, t, t, g);
                break;

            case 4: /* t aliased with g */
                fmpz_poly_set(t, g);
                fmpz_poly_xgcd_modular(r, s, t, f, t);
                break;

            default: FLINT_UNREACHABLE;
        }

        fmpz_poly_mul(tmp, s, f);
        fmpz_poly_mul(res, t, g);
        fmpz_poly_add(res, res, tmp);

        /* Either s * f + t * g == r */
        result = fmpz_poly_equal_fmpz(res, r);

        /* ... Or r is zero with f or g being zero */
        result |= (fmpz_is_zero(r) && (fmpz_poly_is_zero(f) || fmpz_poly_is_zero(g)));
        if (!result)
            TEST_FUNCTION_FAIL(
                    "s * f + t * g != r\n"
                    "type = %d\n"
                    "\n"
                    "s * f + t * g = %{fmpz_poly}\n"
                    "r = %{fmpz}\n"
                    "\n"
                    "f = %{fmpz_poly}\n"
                    "g = %{fmpz_poly}\n"
                    "s = %{fmpz_poly}\n"
                    "t = %{fmpz_poly}\n"
                    "\n"
                    "d = %{fmpz_poly}\n",
                    type, res, r, f, g, s, t, d);

        fmpz_clear(r);
        fmpz_poly_clear(tmp);
        fmpz_poly_clear(res);
        fmpz_poly_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(s);
        fmpz_poly_clear(t);
    }

    TEST_FUNCTION_END(state);
}
