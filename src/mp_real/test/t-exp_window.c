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

/* mp_real_exp_bits on short dyadic arguments of either sign,
   x = +-c 2^-(hi + len) with c of len bits (so |x| in [2^-(hi+1),
   2^-hi)), at precisions beyond the bitwise kernels, where short
   arguments take the bit-burst directly (positive) or 1 / exp(|x|)
   (negative) instead of a reduction that makes them full-length:
   against arb_exp, the result must contain the value and be accurate
   to about n limbs. */

TEST_FUNCTION_START(mp_real_exp_window, state)
{
    slong iter;

    for (iter = 0; iter < 40 * flint_test_multiplier(); iter++)
    {
        slong n, hi, len;
        fmpz_t c;
        mp_real_t x, y;
        arb_t a, r, ref;

        n = 601 + n_randint(state, (iter % 10 == 0) ? 8000 : 1500);
        hi = n_randint(state, 4) ? n_randint(state, 200) : n_randint(state, 3000);
        len = 1 + n_randint(state, n_randint(state, 2) ? 128 : 8192);

        fmpz_init(c);
        mp_real_init(x);
        mp_real_init(y);
        arb_init(a);
        arb_init(r);
        arb_init(ref);

        fmpz_randbits(c, state, len);
        fmpz_abs(c, c);
        fmpz_setbit(c, len - 1);
        mp_real_set_fmpz(x, c);
        mp_real_mul_2exp_si(x, x, -(hi + len));
        if (n_randint(state, 2))
            mp_real_neg(x, x);

        mp_real_get_arb(a, x);
        arb_exp(ref, a, FLINT_BITS * n + 200);
        mp_real_exp_bits(y, x, FLINT_BITS * n);
        mp_real_get_arb(r, y);

        if (!arb_overlaps(r, ref) || arb_rel_accuracy_bits(r) < FLINT_BITS * n - 8)
        {
            flint_printf("FAIL: n = %wd, hi = %wd, len = %wd, sign %d\n", n, hi, len, x->negative ? -1 : 1);
            flint_printf("y = "); arb_printd(r, 30); flint_printf("\n");
            flint_printf("    "); arb_printd(ref, 30); flint_printf("\n");
            flint_printf("accuracy %wd\n", arb_rel_accuracy_bits(r));
            flint_abort();
        }

        fmpz_clear(c);
        mp_real_clear(x);
        mp_real_clear(y);
        arb_clear(a);
        arb_clear(r);
        arb_clear(ref);
    }

    TEST_FUNCTION_END(state);
}
