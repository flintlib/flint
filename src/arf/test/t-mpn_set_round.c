/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_helpers.h"

TEST_FUNCTION_START(_arf_mpn_set_round, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * flint_test_multiplier(); iter++)
    {
        ulong x[20], buf[20], res[20];
        slong xn, prec, rn, exp, bits, i, off;
        fmpz_t X, Q, A;
        int inexact, inexact_ref, in_place;
        arf_rnd_t rnd;

        fmpz_init(X); fmpz_init(Q); fmpz_init(A);

        xn = 1 + n_randint(state, 8);
        _arf_mpn_test_randtest_mantissa(x, xn, state);
        /* leading zero bits, sometimes */
        if (n_randint(state, 2))
            x[xn - 1] >>= n_randint(state, FLINT_BITS);
        if (x[xn - 1] == 0)
            x[xn - 1] = 1;
        /* long runs for the tie and carry cases */
        if (n_randint(state, 3) == 0)
        {
            for (i = 0; i < xn - 1; i++)
                x[i] = n_randint(state, 2) ? ~UWORD(0) : 0;
            if (n_randint(state, 2))
                x[xn - 1] |= ~UWORD(0) >> 1;
        }

        prec = 1 + n_randint(state, 10 * FLINT_BITS);
        rn = (prec + FLINT_BITS - 1) / FLINT_BITS;
        rnd = _arf_mpn_test_rnd_signfree(state);
        fmpz_set_ui_array(X, x, xn);

        /* in place: x inside a buffer, the result at its start (x itself
           when xn = rn, at an offset when xn > rn) */
        in_place = (xn >= rn) && n_randint(state, 2);
        if (in_place)
        {
            off = (xn > rn) ? n_randint(state, 20 - xn + 1) : 0;
            flint_mpn_copyi(buf + off, x, xn);
            inexact = _arf_mpn_set_round(buf, &exp, buf + off, xn, prec, rnd);
            flint_mpn_copyi(res, buf, rn);
        }
        else
            inexact = _arf_mpn_set_round(res, &exp, x, xn, prec, rnd);

        inexact_ref = _arf_mpn_test_round_fmpz(Q, &bits, X, prec, rnd);
        fmpz_set_ui_array(A, res, rn);

        /* value: Q 2^(L + bits - 2 prec) 2^(-xn FLINT_BITS) vs (A / B^rn) 2^exp */
        if (inexact != inexact_ref ||
            !_arf_mpn_test_equal_2exp(Q, fmpz_bits(X) + bits - 2 * prec - xn * FLINT_BITS, A, exp - rn * FLINT_BITS) ||
            !(res[rn - 1] >> (FLINT_BITS - 1)) ||
            (res[0] & ((rn * FLINT_BITS - prec == 0) ? 0 : ((UWORD(1) << (rn * FLINT_BITS - prec)) - 1))) != 0)
        {
            flint_printf("FAIL: xn = %wd, prec = %wd, rnd = %d, in_place = %d\n", xn, prec, (int) rnd, in_place);
            flint_printf("X = "); fmpz_print(X); flint_printf("\n");
            flint_printf("Q = "); fmpz_print(Q); flint_printf(" bits %wd inexact %d\n", bits, inexact_ref);
            flint_printf("A = "); fmpz_print(A); flint_printf(" exp %wd inexact %d\n", exp, inexact);
            flint_abort();
        }

        fmpz_clear(X); fmpz_clear(Q); fmpz_clear(A);
    }

    TEST_FUNCTION_END(state);
}
