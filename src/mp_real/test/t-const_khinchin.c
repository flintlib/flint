/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "mp_real.h"

/* Khinchin's constant with forced parameters (the logarithm rows below
   K0, the Bernoulli columns below n1, the rest summed directly, which
   covers every combination of the three parts, all-Bernoulli and
   almost all-direct), on 1 to 4 threads with 1 to 4 Bernoulli chunks
   at once (each with the threads left over inside, from 1500 limbs
   now and then, where the iterator uses them too): the balls must overlap a
   reference at a higher precision and other parameters, contain the
   160-digit reference (computed with mpmath) as far as it reaches, and
   be accurate to within a few bits of n limbs. */

static const char * khinchin_ref = "2.685452001065306445309714835481795693820382293994462953051152345557218859537152002801141174931847697995153465905288090082897677716410963051792533483259668381852";

TEST_FUNCTION_START(mp_real_const_khinchin, state)
{
    slong iter;
    arb_t ref160;

    arb_init(ref160);
    arb_set_str(ref160, khinchin_ref, 560);
    arb_add_error_2exp_si(ref160, -525);

    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong n, K0, n1, nb, acc;
        mp_real_t x, y;
        arb_t a, b;

        if (iter % 25 == 24)
            n = 1500 + n_randint(state, 300);
        else
            n = 1 + n_randint(state, (iter % 10 == 0) ? 200 : 40);

        if (n_randint(state, 4) == 0)
            K0 = 0;
        else
            K0 = 3 + n_randint(state, (n_randint(state, 4) == 0) ? 60 : 12);

        switch (n_randint(state, 5))
        {
            case 0:
                n1 = 0;                     /* default */
                break;
            case 1:
                n1 = WORD_MAX;              /* all Bernoulli */
                break;
            case 2:
                n1 = 1;                     /* as few Bernoulli as allowed */
                break;
            default:
                n1 = 1 + n_randint(state, FLINT_BITS * n);
        }

        flint_set_num_threads(1 + n_randint(state, 4));
        nb = n_randint(state, 5);

        mp_real_init(x);
        mp_real_init(y);
        arb_init(a);
        arb_init(b);

        _mp_real_const_khinchin_tune(x, n, K0, (n >= 1500) ? 0 : n1, nb);
        mp_real_get_arb(a, x);

        _mp_real_const_khinchin_tune(y, n + 1 + n_randint(state, 3),
            0, 0, 0);
        mp_real_get_arb(b, y);

        if (!arb_overlaps(a, b) || !arb_overlaps(a, ref160))
        {
            flint_printf("FAIL: overlap, n = %wd, K0 = %wd, n1 = %wd, nb = %wd\n", n, K0, n1, nb);
            flint_printf("a = "); arb_printn(a, 50, 0); flint_printf("\n");
            flint_printf("b = "); arb_printn(b, 50, 0); flint_printf("\n");
            flint_abort();
        }

        acc = arb_rel_accuracy_bits(a);
        if (acc < FLINT_BITS * n - 40)
        {
            flint_printf("FAIL: accuracy %wd, n = %wd, K0 = %wd, n1 = %wd, nb = %wd\n",
                acc, n, K0, n1, nb);
            flint_abort();
        }

        mp_real_clear(x);
        mp_real_clear(y);
        arb_clear(a);
        arb_clear(b);
    }

    arb_clear(ref160);
    flint_set_num_threads(1);

    TEST_FUNCTION_END(state);
}
