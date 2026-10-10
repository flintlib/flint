/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"
#include "arb_poly.h"

static const ulong known_values[] =
{
    UWORD(2147483629),
    UWORD(1073742093),
    UWORD(1342248677),
    UWORD(3319936736),
    UWORD(2947821228),
    UWORD(1019513834),
    UWORD(3324951530),
    UWORD(1995039408),
    UWORD(3505683295),
    UWORD(3567639420),
    UWORD(394942914)
};

/* The previous implementation: ball arithmetic with rounding. */
static void
swinnerton_dyer_arb(fmpz_poly_t poly, ulong n)
{
    slong N = (WORD(1) << n);

    fmpz_poly_fit_length(poly, N + 1);

    if (n == 0)
    {
        fmpz_zero(poly->coeffs);
        fmpz_one(poly->coeffs + 1);
    }
    else
    {
        arb_poly_t t;
        arb_poly_init(t);
        arb_poly_swinnerton_dyer_ui(t, n, 0);
        if (!_arb_vec_get_unique_fmpz_vec(poly->coeffs, t->coeffs, t->length))
            flint_throw(FLINT_ERROR, "arb reference failed for n = %wu\n", n);
        arb_poly_clear(t);
    }

    _fmpz_poly_set_length(poly, N + 1);
}

TEST_FUNCTION_START(fmpz_poly_swinnerton_dyer, state)
{
    fmpz_poly_t S, T, U;
    ulong r;
    slong n, i;

    fmpz_poly_init(S);
    fmpz_poly_init(T);
    fmpz_poly_init(U);

    /* fixed small cases */
    fmpz_poly_swinnerton_dyer(S, 0);
    fmpz_poly_set_str(T, "2  0 1");
    FLINT_TEST(fmpz_poly_equal(S, T));

    fmpz_poly_swinnerton_dyer(S, 1);
    fmpz_poly_set_str(T, "3  -2 0 1");
    FLINT_TEST(fmpz_poly_equal(S, T));

    fmpz_poly_swinnerton_dyer(S, 2);
    fmpz_poly_set_str(T, "5  1 0 -10 0 1");
    FLINT_TEST(fmpz_poly_equal(S, T));

    fmpz_poly_swinnerton_dyer(S, 3);
    fmpz_poly_set_str(T, "9  576 0 -960 0 352 0 -40 0 1");
    FLINT_TEST(fmpz_poly_equal(S, T));

    /* known evaluations modulo a prime, independent of the arb code */
    for (n = 0; n <= 10; n++)
    {
        fmpz_poly_swinnerton_dyer(S, n);
        r = fmpz_poly_evaluate_mod(S, UWORD(2147483629), UWORD(4294967291));

        if (r != known_values[n])
        {
            flint_printf("ERROR: wrong evaluation of S_%wd\n", n);
            fflush(stdout);
            flint_abort();
        }
    }

    /* exact agreement with the arb-based implementation; the larger n
       exercise the divide-and-conquer branch of the Taylor shift */
    for (n = 0; n <= 11; n++)
    {
        fmpz_poly_swinnerton_dyer(S, n);
        swinnerton_dyer_arb(T, n);

        if (!fmpz_poly_equal(S, T))
        {
            flint_printf("FAIL: new and arb implementations differ, n = %wd\n", n);
            flint_printf("new: "); fmpz_poly_print_pretty(S, "x"); flint_printf("\n");
            flint_printf("arb: "); fmpz_poly_print_pretty(T, "x"); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }
    }

    /* structure: monic, degree 2^n, even */
    for (n = 1; n <= 11; n++)
    {
        fmpz_poly_swinnerton_dyer(S, n);

        FLINT_TEST(fmpz_poly_degree(S) == (WORD(1) << n));
        FLINT_TEST(fmpz_is_one(fmpz_poly_lead(S)));

        for (i = 1; i < S->length; i += 2)
            FLINT_TEST(fmpz_is_zero(S->coeffs + i));
    }

    /* the output polynomial is reused with decreasing and increasing n */
    for (n = 9; n >= 0; n--)
    {
        fmpz_poly_swinnerton_dyer(S, n);
        swinnerton_dyer_arb(T, n);
        FLINT_TEST(fmpz_poly_equal(S, T));
    }

    for (n = 0; n <= 9; n++)
    {
        fmpz_poly_swinnerton_dyer(S, n);
        swinnerton_dyer_arb(T, n);
        FLINT_TEST(fmpz_poly_equal(S, T));
    }

    /* the vector and polynomial interfaces agree */
    for (n = 0; n <= 8; n++)
    {
        fmpz_poly_swinnerton_dyer(S, n);

        fmpz_poly_fit_length(U, (WORD(1) << n) + 1);
        _fmpz_poly_swinnerton_dyer(U->coeffs, n);
        _fmpz_poly_set_length(U, (WORD(1) << n) + 1);

        FLINT_TEST(fmpz_poly_equal(S, U));
    }

    fmpz_poly_clear(S);
    fmpz_poly_clear(T);
    fmpz_poly_clear(U);

    TEST_FUNCTION_END(state);
}
