/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "radix.h"
#include "gr.h"

TEST_FUNCTION_START(radix_integer, state)
{
    slong iter;

    /* get/set_limb */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        radix_integer_t x, y, z;
        ulong c, d, e;
        slong i;

        radix_init_randtest(radix, state);
        radix_integer_init(x, radix);
        radix_integer_init(y, radix);
        radix_integer_init(z, radix);

        radix_integer_randtest_limbs(x, state, 5, radix);
        radix_integer_randtest_limbs(y, state, 5, radix);

        radix_randtest_limbs(&c, state, 1, radix);
        i = n_randint(state, 5);

        radix_integer_set_limb(y, x, i, c, radix);
        if (n_randint(state, 2))
            radix_integer_neg(y, y, radix);

        radix_integer_set(z, x, radix);
        if (n_randint(state, 2))
            radix_integer_neg(z, z, radix);
        radix_integer_set_limb(z, z, i, c, radix);
        if (n_randint(state, 2))
            radix_integer_neg(z, z, radix);

        d = radix_integer_get_limb(y, i, radix);
        e = radix_integer_get_limb(z, i, radix);

        if (d != c || e != c || !radix_integer_is_normalised(y, radix) ||
                                !radix_integer_is_normalised(z, radix))
        {
            flint_printf("FAIL: get/set_limb\n");
            flint_printf("x = (%wd %wd %{ulong*})\n", x->size, x->alloc, x->d, radix_integer_size(x, radix));
            flint_printf("y = (%wd %wd %{ulong*})\n", y->size, y->alloc, y->d, radix_integer_size(y, radix));
            flint_printf("z = (%wd %wd %{ulong*})\n", z->size, z->alloc, z->d, radix_integer_size(z, radix));
            flint_printf("i = %wd\n", i);
            flint_printf("c = %wu\n", c);
            flint_printf("d = %wu\n", d);
            flint_printf("e = %wu\n", e);
            flint_abort();
        }

        radix_integer_clear(x, radix);
        radix_integer_clear(y, radix);
        radix_integer_clear(z, radix);
        radix_clear(radix);
    }

    /* lshift/rshift_limbs */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        radix_integer_t x, y, z;
        slong i, n;

        radix_init_randtest(radix, state);
        radix_integer_init(x, radix);
        radix_integer_init(y, radix);
        radix_integer_init(z, radix);

        radix_integer_randtest_limbs(x, state, 5, radix);
        radix_integer_randtest_limbs(y, state, 5, radix);

        if (n_randint(state, 2))
        {
            n = n_randint(state, 5);

            if (n_randint(state, 2))
            {
                radix_integer_set(y, x, radix);
                radix_integer_lshift_limbs(y, y, n, radix);
            }
            else
            {
                radix_integer_lshift_limbs(y, x, n, radix);
            }

            for (i = 0; i < n; i++)
                radix_integer_set_limb(y, y, i, 1, radix);

            if (n_randint(state, 2))
            {
                radix_integer_set(z, y, radix);
                radix_integer_rshift_limbs(z, z, n, radix);
            }
            else
            {
                radix_integer_rshift_limbs(z, y, n, radix);
            }
        }
        else
        {
            n = radix_integer_valuation_limbs(x, radix);
            radix_integer_rshift_limbs(y, x, n, radix);
            radix_integer_lshift_limbs(z, y, n, radix);
        }

        if (!radix_integer_equal(x, z, radix))
        {
            flint_printf("FAIL: lshift/rshift\n");
            flint_printf("x = (%wd %wd %{ulong*})\n", x->size, x->alloc, x->d, radix_integer_size(x, radix));
            flint_printf("y = (%wd %wd %{ulong*})\n", y->size, y->alloc, y->d, radix_integer_size(y, radix));
            flint_printf("z = (%wd %wd %{ulong*})\n", z->size, z->alloc, z->d, radix_integer_size(z, radix));
            flint_printf("n = %wd\n", n);
            flint_abort();
        }

        radix_integer_clear(x, radix);
        radix_integer_clear(y, radix);
        radix_integer_clear(z, radix);
        radix_clear(radix);
    }

    /* mod_limbs */
    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        radix_integer_t x, y, z;
        slong n;

        radix_init_randtest(radix, state);
        radix_integer_init(x, radix);
        radix_integer_init(y, radix);
        radix_integer_init(z, radix);

        radix_integer_randtest_limbs(x, state, 5, radix);
        radix_integer_randtest_limbs(y, state, 5, radix);
        radix_integer_randtest_limbs(z, state, 5, radix);

        n = n_randint(state, 5);

        radix_integer_mod_limbs(y, x, n, radix);
        radix_integer_lshift_limbs(z, z, n, radix);
        radix_integer_add(z, z, x, radix);
        radix_integer_mod_limbs(z, z, n, radix);

        if (!radix_integer_equal(y, z, radix) ||
            !radix_integer_is_normalised(y, radix) ||
            !radix_integer_is_normalised(z, radix) ||
            radix_integer_sgn(y, radix) < 0)
        {
            flint_printf("FAIL: mod_limbs\n");
            flint_printf("radix %wu ^ %u = %wu\n", DIGIT_RADIX(radix), radix->exp, LIMB_RADIX(radix));
            flint_printf("x = (%wd %wd %{ulong*})\n", x->size, x->alloc, x->d, radix_integer_size(x, radix));
            flint_printf("y = (%wd %wd %{ulong*})\n", y->size, y->alloc, y->d, radix_integer_size(y, radix));
            flint_printf("z = (%wd %wd %{ulong*})\n", z->size, z->alloc, z->d, radix_integer_size(z, radix));
            flint_printf("n = %wd\n", n);
            flint_abort();
        }

        radix_integer_clear(x, radix);
        radix_integer_clear(y, radix);
        radix_integer_clear(z, radix);
        radix_clear(radix);
    }

    /* smod_limbs */
    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        radix_integer_t x, y, z;
        slong n;

        radix_init_randtest(radix, state);
        radix_integer_init(x, radix);
        radix_integer_init(y, radix);
        radix_integer_init(z, radix);

        radix_integer_randtest_limbs(x, state, 5, radix);
        radix_integer_randtest_limbs(y, state, 5, radix);
        radix_integer_randtest_limbs(z, state, 5, radix);

        n = n_randint(state, 5);

        radix_integer_smod_limbs(y, x, n, radix);
        radix_integer_lshift_limbs(z, z, n, radix);
        radix_integer_add(z, z, x, radix);
        radix_integer_smod_limbs(z, z, n, radix);

        if (!radix_integer_equal(y, z, radix) ||
            !radix_integer_is_normalised(y, radix) ||
            !radix_integer_is_normalised(z, radix) ||
            radix_integer_cmpabs(y, x, radix) > 0)
        {
            flint_printf("FAIL: smod_limbs\n");
            flint_printf("radix %wu ^ %u = %wu\n", DIGIT_RADIX(radix), radix->exp, LIMB_RADIX(radix));
            flint_printf("x = (%wd %wd %{ulong*})\n", x->size, x->alloc, x->d, radix_integer_size(x, radix));
            flint_printf("y = (%wd %wd %{ulong*})\n", y->size, y->alloc, y->d, radix_integer_size(y, radix));
            flint_printf("z = (%wd %wd %{ulong*})\n", z->size, z->alloc, z->d, radix_integer_size(z, radix));
            flint_printf("n = %wd\n", n);
            flint_abort();
        }

        radix_integer_clear(x, radix);
        radix_integer_clear(y, radix);
        radix_integer_clear(z, radix);
        radix_clear(radix);
    }

    /* mullow */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        radix_integer_t x, y, z, w;
        slong n, zrn;

        radix_init_randtest(radix, state);
        radix_integer_init(x, radix);
        radix_integer_init(y, radix);
        radix_integer_init(z, radix);
        radix_integer_init(w, radix);

        radix_integer_randtest_limbs(x, state, 4, radix);
        radix_integer_randtest_limbs(y, state, 4, radix);
        radix_integer_randtest_limbs(z, state, 4, radix);

        n = n_randint(state, 4);

        switch (n_randint(state, 5))
        {
            case 0:
                radix_integer_mullow_limbs(z, x, y, n, radix);
                break;
            case 1:
                radix_integer_set(y, x, radix);
                radix_integer_mullow_limbs(z, x, x, n, radix);
                break;
            case 2:
                radix_integer_set(z, x, radix);
                radix_integer_mullow_limbs(z, z, y, n, radix);
                break;
            case 3:
                radix_integer_set(z, y, radix);
                radix_integer_mullow_limbs(z, x, z, n, radix);
                break;
            default:
                radix_integer_set(y, x, radix);
                radix_integer_set(z, x, radix);
                radix_integer_mullow_limbs(z, z, z, n, radix);
                break;
        }

        zrn = radix_integer_size(z, radix);
        radix_integer_smod_limbs(z, z, n, radix);

        radix_integer_mul(w, x, y, radix);
        radix_integer_smod_limbs(w, w, n, radix);

        if (zrn > n || !radix_integer_equal(w, z, radix) ||
            !radix_integer_is_normalised(w, radix))
        {
            flint_printf("FAIL: mullow_limbs\n");
            flint_printf("radix %wu ^ %u = %wu\n", DIGIT_RADIX(radix), radix->exp, LIMB_RADIX(radix));
            flint_printf("x = (%wd %wd %{ulong*})\n", x->size, x->alloc, x->d, radix_integer_size(x, radix));
            flint_printf("y = (%wd %wd %{ulong*})\n", y->size, y->alloc, y->d, radix_integer_size(y, radix));
            flint_printf("z = (%wd %wd %{ulong*})\n", z->size, z->alloc, z->d, radix_integer_size(z, radix));
            flint_printf("w = (%wd %wd %{ulong*})\n", w->size, w->alloc, w->d, radix_integer_size(w, radix));
            flint_printf("n = %wd\n", n);
            flint_abort();
        }

        radix_integer_clear(x, radix);
        radix_integer_clear(y, radix);
        radix_integer_clear(z, radix);
        radix_integer_clear(w, radix);
        radix_clear(radix);
    }

    /* invmod */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        radix_integer_t x, y, z;
        slong n;
        int invertible;

        radix_init_randtest(radix, state);
        radix_integer_init(x, radix);
        radix_integer_init(y, radix);
        radix_integer_init(z, radix);

        radix_integer_randtest_limbs(x, state, 4, radix);
        radix_integer_randtest_limbs(y, state, 4, radix);
        radix_integer_randtest_limbs(z, state, 4, radix);

        n = n_randint(state, 4);

        if (n_randint(state, 2))
        {
            invertible = radix_integer_invmod_limbs(y, x, n, radix);
        }
        else
        {
            radix_integer_set(y, x, radix);
            invertible = radix_integer_invmod_limbs(y, y, n, radix);
        }

        if (invertible)
        {
            radix_integer_mullow_limbs(z, x, y, n, radix);

            if ((n == 0 && !radix_integer_is_zero(y, radix)) ||
                (n != 0 && (
                    !radix_integer_is_one(z, radix) ||
                    !radix_integer_is_normalised(y, radix) ||
                    radix_integer_sgn(y, radix) != radix_integer_sgn(x, radix))))
            {
                flint_printf("FAIL: invmod_limbs\n");
                flint_printf("radix %wu ^ %u = %wu\n", DIGIT_RADIX(radix), radix->exp, LIMB_RADIX(radix));
                flint_printf("x = (%wd %wd %{ulong*})\n", x->size, x->alloc, x->d, radix_integer_size(x, radix));
                flint_printf("y = (%wd %wd %{ulong*})\n", y->size, y->alloc, y->d, radix_integer_size(y, radix));
                flint_printf("z = (%wd %wd %{ulong*})\n", z->size, z->alloc, z->d, radix_integer_size(z, radix));
                flint_printf("n = %wd\n", n);
                flint_abort();
            }
        }

        radix_integer_clear(x, radix);
        radix_integer_clear(y, radix);
        radix_integer_clear(z, radix);
        radix_clear(radix);
    }

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        gr_ctx_t ctx;

        radix_init_randtest(radix, state);
        gr_ctx_init_radix_integer(ctx, DIGIT_RADIX(radix), radix->exp);

        gr_test_ring(ctx, 10, 0);

        gr_ctx_clear(ctx);
        radix_clear(radix);
    }

    TEST_FUNCTION_END(state);
}
