/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_helpers.h"

/* an approximation (a, n) with guard g of an exact value (the same with
   limbs below, perturbed by less than 2^e ulps of g in either direction,
   or in the known direction of side) rounds like the exact value
   whenever the window certifies it */
TEST_FUNCTION_START(_arf_mpn_round_guarded, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * flint_test_multiplier(); iter++)
    {
        ulong a[12], full[14], res[12];
        slong n, prec, rn, exp, i, d;
        fmpz_t A, P, Q, R;
        int e, r, inexact_ref, sign, side;
        arf_rnd_t rnd;

        fmpz_init(A); fmpz_init(P); fmpz_init(Q); fmpz_init(R);

        n = 1 + n_randint(state, 8);
        _arf_mpn_test_randtest_mantissa(a, n, state);
        if (n_randint(state, 3) == 0)
            for (i = 0; i < n - 1; i++)
                a[i] = n_randint(state, 2) ? ~UWORD(0) : 0;
        prec = 1 + n_randint(state, n * FLINT_BITS);
        rn = (prec + FLINT_BITS - 1) / FLINT_BITS;
        rnd = _arf_mpn_test_rnd_signfree(state);
        e = n_randint(state, 20);
        side = (int) n_randint(state, 3) - 1;
        /* an approximation exact at prec bits, now and then */
        if (n_randint(state, 4) == 0 && rn * FLINT_BITS - prec != 0)
            a[n - rn] &= ~UWORD(0) << (rn * FLINT_BITS - prec);

        /* the approximation A = (a, g) as an integer of n + 1 limbs, and
           the exact value P = A B + low with |P - A B| < 2^e B (two
           limbs below a: the guard g and one more) */
        full[0] = n_randtest(state);
        full[1] = n_randint(state, 2) ? n_randtest(state) : (n_randint(state, 2) ? 0 : ~UWORD(0));
        flint_mpn_copyi(full + 2, a, n);
        fmpz_set_ui_array(A, full + 1, n + 1);
        fmpz_set_ui_array(P, full, n + 2);
        /* perturb by d with |d| < 2^e ulps of the guard limb, i.e. |d B| */
        if (side != 0)
        {
            /* P = A B + R with 0 < |R| < 2^e B, strictly in the direction
               of side */
            fmpz_mul_2exp(P, A, FLINT_BITS);
            fmpz_randtest_unsigned(R, state, e + FLINT_BITS);
            if (fmpz_is_zero(R))
                fmpz_one(R);
            if (side < 0)
                fmpz_sub(P, P, R);
            else
                fmpz_add(P, P, R);
        }
        else
        {
            d = (e == 0) ? 0 : (slong) n_randint(state, UWORD(1) << e);
            sign = n_randint(state, 2);
            fmpz_set_ui(R, d);
            fmpz_mul_2exp(R, R, FLINT_BITS);
            if (sign)
                fmpz_sub(P, P, R);
            else
                fmpz_add(P, P, R);
        }
        if (fmpz_sgn(P) <= 0 || fmpz_bits(P) > (n + 2) * FLINT_BITS || fmpz_bits(P) < (n + 2) * FLINT_BITS - 1)
        {
            /* the perturbation crossed a power of two upwards (or fell
               far): skip; falling just below one is the case where the
               approximation is a power of two, which is handled */
            fmpz_clear(A); fmpz_clear(P); fmpz_clear(Q); fmpz_clear(R);
            continue;
        }

        r = _arf_mpn_round_guarded(res, &exp, a, n, full[1], e, prec, rnd, side);

        if (r >= 0)
        {
            slong bits;
            fmpz_t Rr;
            fmpz_init(Rr);

            inexact_ref = _arf_mpn_test_round_fmpz(Q, &bits, P, prec, rnd);
            fmpz_set_ui_array(Rr, res, rn);

            /* value: Q 2^(L + bits - 2 prec) with L the bit length of P,
               vs res 2^(exp + (n + 2 - rn) FLINT_BITS) (the approximation
               has (n + 2) FLINT_BITS bits) */
            if (inexact_ref != 1 || exp < -1 || exp > 1 || !(res[rn - 1] >> (FLINT_BITS - 1)) ||
                !_arf_mpn_test_equal_2exp(Q, fmpz_bits(P) + bits - 2 * prec, Rr, exp + (n + 2 - rn) * FLINT_BITS))
            {
                flint_printf("FAIL: n = %wd, prec = %wd, rnd = %d, e = %d, side = %d\n", n, prec, (int) rnd, e, side);
                flint_printf("P = "); fmpz_print(P); flint_printf("\n");
                flint_printf("Q = "); fmpz_print(Q); flint_printf(" bits %wd inexact %d\n", bits, inexact_ref);
                flint_printf("res = "); fmpz_print(Rr); flint_printf(" exp %wd\n", exp);
                flint_abort();
            }
            fmpz_clear(Rr);
        }

        fmpz_clear(A); fmpz_clear(P); fmpz_clear(Q); fmpz_clear(R);
    }

    TEST_FUNCTION_END(state);
}
