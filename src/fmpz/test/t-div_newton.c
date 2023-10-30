/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
test_div_q(void (*my_f)(fmpz_t, const fmpz_t, const fmpz_t),
           void (*reference_f)(fmpz_t, const fmpz_t, const fmpz_t), flint_rand_t state, const char * descr)
{
    fmpz_t a, b, q, q2;
    int aliasing;
    slong bits;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(q);
    fmpz_init(q2);

    if (n_randint(state, 100) == 0)
        bits = 100000;
    else
        bits = 200;

    if (reference_f == fmpz_divexact)
    {
        fmpz_randtest(a, state, bits);
        fmpz_randtest_not_zero(b, state, bits);
        fmpz_mul(a, a, b);
    }
    else
    {
        fmpz_randtest(a, state, bits);
        fmpz_randtest_not_zero(b, state, bits);

        if (n_randint(state, 4) == 0)
        {
            fmpz_mul(a, a, b);
            fmpz_add_si(a, a, (slong) n_randint(state, 5) - 2);
        }
    }

    aliasing = n_randint(state, 3);

    if (aliasing == 0)
    {
        my_f(q, a, b);
    }
    else if (aliasing == 1)
    {
        fmpz_set(q, a);
        my_f(q, q, b);
    }
    else
    {
        fmpz_set(q, b);
        my_f(q, a, q);
    }

    reference_f(q2, a, b);

    if (!fmpz_equal(q, q2) || !_fmpz_is_canonical(q))
    {
        flint_printf("FAIL: %s\n", descr);
        flint_printf("aliasing = %d\n", aliasing);
        flint_printf("a = "); fmpz_print(a); flint_printf("\n");
        flint_printf("b = "); fmpz_print(b); flint_printf("\n");
        flint_printf("q = "); fmpz_print(q); flint_printf("\n");
        flint_printf("q2 = "); fmpz_print(q2); flint_printf("\n");
        fflush(stdout);
        flint_abort();
    }

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(q);
    fmpz_clear(q2);
}

void
test_div_qr(void (*my_f)(fmpz_t, fmpz_t, const fmpz_t, const fmpz_t),
           void (*reference_f)(fmpz_t, fmpz_t, const fmpz_t, const fmpz_t), flint_rand_t state, const char * descr)
{
    fmpz_t a, b, q, r, q2, r2;
    int aliasing;
    slong bits;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(q);
    fmpz_init(r);
    fmpz_init(q2);
    fmpz_init(r2);

    if (n_randint(state, 100) == 0)
        bits = 100000;
    else
        bits = 200;

    fmpz_randtest(a, state, bits);
    fmpz_randtest_not_zero(b, state, bits);

    if (n_randint(state, 4) == 0)
    {
        fmpz_mul(a, a, b);
        fmpz_add_si(a, a, (slong) n_randint(state, 5) - 2);
    }

    aliasing = n_randint(state, 5);

    if (aliasing == 0)
    {
        my_f(q, r, a, b);
    }
    else if (aliasing == 1)
    {
        fmpz_set(q, a);
        my_f(q, r, q, b);
    }
    else if (aliasing == 2)
    {
        fmpz_set(q, b);
        my_f(q, r, a, q);
    }
    else if (aliasing == 3)
    {
        fmpz_set(r, a);
        my_f(q, r, r, b);
    }
    else
    {
        fmpz_set(r, b);
        my_f(q, r, a, r);
    }

    reference_f(q2, r2, a, b);

    if (!fmpz_equal(q, q2) || !fmpz_equal(r, r2) || !_fmpz_is_canonical(q) || !_fmpz_is_canonical(r))
    {
        flint_printf("FAIL: %s\n", descr);
        flint_printf("aliasing = %d\n", aliasing);
        flint_printf("a = "); fmpz_print(a); flint_printf("\n");
        flint_printf("b = "); fmpz_print(b); flint_printf("\n");
        flint_printf("q = "); fmpz_print(q); flint_printf("\n");
        flint_printf("r = "); fmpz_print(r); flint_printf("\n");
        flint_printf("q2 = "); fmpz_print(q2); flint_printf("\n");
        flint_printf("r2 = "); fmpz_print(r2); flint_printf("\n");
        fflush(stdout);
        flint_abort();
    }

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(q);
    fmpz_clear(r);
    fmpz_clear(q2);
    fmpz_clear(r2);
}

/* these currently don't exist in flint */
void
fmpz_tdiv_r(fmpz_t r, const fmpz_t a, const fmpz_t b)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_tdiv_qr(t, r, a, b);
    fmpz_clear(t);
}

void
fmpz_cdiv_r(fmpz_t r, const fmpz_t a, const fmpz_t b)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_cdiv_qr(t, r, a, b);
    fmpz_clear(t);
}

TEST_FUNCTION_START(fmpz_div_newton, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        test_div_q(_fmpz_tdiv_q_newton, fmpz_tdiv_q, state, "tdiv_q");
        test_div_q(_fmpz_fdiv_q_newton, fmpz_fdiv_q, state, "fdiv_q");
        test_div_q(_fmpz_cdiv_q_newton, fmpz_cdiv_q, state, "cdiv_q");

        test_div_q(_fmpz_tdiv_r_newton, fmpz_tdiv_r, state, "tdiv_r");
        test_div_q(_fmpz_fdiv_r_newton, fmpz_fdiv_r, state, "fdiv_r");
        test_div_q(_fmpz_cdiv_r_newton, fmpz_cdiv_r, state, "cdiv_r");

        test_div_q(_fmpz_mod_newton, fmpz_mod, state, "mod");

        test_div_qr(_fmpz_tdiv_qr_newton, fmpz_tdiv_qr, state, "tdiv_qr");
        test_div_qr(_fmpz_fdiv_qr_newton, fmpz_fdiv_qr, state, "fdiv_qr");
        test_div_qr(_fmpz_cdiv_qr_newton, fmpz_cdiv_qr, state, "cdiv_qr");

        test_div_q(_fmpz_divexact_newton, fmpz_divexact, state, "divexact");
    }

    TEST_FUNCTION_END(state);
}
