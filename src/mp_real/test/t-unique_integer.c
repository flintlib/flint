/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "arb.h"
#include "mp_real.h"
#include "mp_real/impl.h"

/* mp_real_unique_integer against arb_get_unique_fmpz (exact here: radii
   below 2^30 ulps convert to arb without rounding), on balls built
   around integers with fractions and radii at the edges (fractions of
   all-zero or all-ones limbs, radius equal to the fraction or its
   complement), aliased or not; the result is exact and normalized, and
   res is unchanged on failure. */

static ulong
_ui_special(flint_rand_t state)
{
    switch (n_randint(state, 6))
    {
        case 0: return 0;
        case 1: return ~UWORD(0);
        case 2: return n_randint(state, 4);
        case 3: return ~UWORD(0) - n_randint(state, 4);
        case 4: return UWORD(1) << (FLINT_BITS - 1);
        default: return n_randtest(state);
    }
}

TEST_FUNCTION_START(mp_real_unique_integer, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * flint_test_multiplier(); iter++)
    {
        mp_real_t x, y, s;
        arb_t a, b;
        fmpz_t f, g;
        slong size, f_limbs, i;
        int ok1, ok2, alias = n_randint(state, 2);

        mp_real_init(x); mp_real_init(y); mp_real_init(s);
        arb_init(a); arb_init(b);
        fmpz_init(f); fmpz_init(g);

        /* x = (d, size) B^(exp - size), exp - size = -f_limbs (or >= 0) */
        size = n_randint(state, 6);
        mp_real_fit_length(x, size + 1);
        for (i = 0; i < size; i++)
            x->d[i] = _ui_special(state);
        if (size > 0 && x->d[size - 1] == 0)
            x->d[size - 1] = 1 + n_randint(state, 3);
        x->size = size;
        f_limbs = (slong) n_randint(state, 8) - 1;      /* -1 .. 6 */
        x->exp = size - f_limbs;
        x->negative = (size > 0) && n_randint(state, 2);
        switch (n_randint(state, 5))
        {
            case 0: x->err = 0; break;
            case 1: x->err = (size > 0) ? x->d[0] : 1; break;          /* F = err */
            case 2: x->err = (size > 0) ? -x->d[0] : 1; break;         /* F + err = B */
            case 3: x->err = 1 + n_randint(state, 3); break;
            default: x->err = n_randint(state, UWORD(1) << 29); break;
        }
        x->err &= (UWORD(1) << 29) - 1;    /* exact in arb (mag) */
        if (x->err == 0)
            _mp_real_norm(x);              /* exact values: no low zero limbs */

        mp_real_get_arb(a, x);
        ok1 = arb_get_unique_fmpz(f, a);

        mp_real_set_si(y, -7);             /* to see that a failure leaves y */
        if (alias)
        {
            mp_real_set(s, x);
            ok2 = mp_real_unique_integer(s, s);
            if (ok2)
                mp_real_swap(y, s);
            else
            {
                mp_real_get_arb(b, s);
                if (!arb_equal(a, b) || s->size != x->size || s->exp != x->exp)
                {
                    flint_printf("FAIL: aliased input changed on failure\n");
                    flint_abort();
                }
            }
        }
        else
            ok2 = mp_real_unique_integer(y, x);

        if (ok1 != ok2)
        {
            flint_printf("FAIL: success (arb %d, mp_real %d), alias %d\n", ok1, ok2, alias);
            mp_real_print(x);
            flint_abort();
        }
        if (ok2)
        {
            mp_real_get_arb(b, y);
            if (!arb_is_exact(b) || !arb_get_unique_fmpz(g, b) || !fmpz_equal(f, g)
                || (y->size > 0 && (y->d[0] == 0 || y->d[y->size - 1] == 0))
                || (y->size == 0 && y->negative))
            {
                flint_printf("FAIL: value, alias %d\n", alias);
                mp_real_print(x); mp_real_print(y);
                fmpz_print(f); flint_printf("\n");
                flint_abort();
            }
        }
        else if (!alias && !(y->size == 1 && y->negative && y->d[0] == 7 && y->err == 0))
        {
            flint_printf("FAIL: res changed on failure\n");
            flint_abort();
        }

        mp_real_clear(x); mp_real_clear(y); mp_real_clear(s);
        arb_clear(a); arb_clear(b);
        fmpz_clear(f); fmpz_clear(g);
    }

    /* random balls of all sizes */
    for (iter = 0; iter < 20000 * flint_test_multiplier(); iter++)
    {
        mp_real_t x, y;
        arb_t a, b;
        fmpz_t f, g;
        int ok1, ok2;

        mp_real_init(x); mp_real_init(y);
        arb_init(a); arb_init(b);
        fmpz_init(f); fmpz_init(g);

        mp_real_randtest(x, state, 1 + n_randint(state, 10), n_randint(state, 2));
        x->exp += (slong) n_randint(state, 12) - 6;
        x->err &= (UWORD(1) << 29) - 1;
        if (x->err == 0)
            _mp_real_norm(x);
        mp_real_get_arb(a, x);
        ok1 = arb_get_unique_fmpz(f, a);
        ok2 = mp_real_unique_integer(y, x);
        if (ok1 != ok2 || (ok2 && (mp_real_get_arb(b, y), !arb_get_unique_fmpz(g, b)
                || !fmpz_equal(f, g))))
        {
            flint_printf("FAIL: random ball\n");
            mp_real_print(x);
            flint_abort();
        }

        mp_real_clear(x); mp_real_clear(y);
        arb_clear(a); arb_clear(b);
        fmpz_clear(f); fmpz_clear(g);
    }

    TEST_FUNCTION_END(state);
}
