/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_next_calkin_wilf, state)
{
    slong i;
    fmpq_t r, ans;


    fmpq_init(r);
    fmpq_init(ans);

    fmpq_set_si(r, 0, 1);
    fmpq_set_si(ans, 43, 205);
    for (i = 0; i < 10000; i++)
        fmpq_next_calkin_wilf(r, r);

    if (!fmpq_equal(r, ans))
    {
        flint_printf("FAIL: enum from 0\n");
        fmpq_print(r);
        flint_printf("\n");
        fmpq_print(ans);
        flint_printf("\n");
        fflush(stdout);
        flint_abort();
    }

    fmpq_set_si(r, 0, 1);
    fmpq_set_si(ans, -43, 162);
    for (i = 0; i < 10000; i++)
        fmpq_next_signed_calkin_wilf(r, r);

    if (!fmpq_equal(r, ans))
    {
        flint_printf("FAIL: signed enum from 0\n");
        fmpq_print(r);
        flint_printf("\n");
        fmpq_print(ans);
        flint_printf("\n");
        fflush(stdout);
        flint_abort();
    }

    fmpz_set_str(fmpq_numref(r), "18446744073709551615", 10);
    fmpz_set_str(fmpq_denref(r), "18446744073709551616", 10);
    fmpz_set_str(fmpq_numref(ans), "498062089990157893394", 10);
    fmpz_set_str(fmpq_denref(ans), "700976274800962961073", 10);
    for (i = 0; i < 1000; i++)
        fmpq_next_calkin_wilf(r, r);

    if (!fmpq_equal(r, ans))
    {
        flint_printf("FAIL: enum from 2^64\n");
        fmpq_print(r);
        flint_printf("\n");
        fmpq_print(ans);
        flint_printf("\n");
        fflush(stdout);
        flint_abort();
    }

    fmpz_set_str(fmpq_numref(r), "18446744073709551615", 10);
    fmpz_set_str(fmpq_denref(r), "18446744073709551616", 10);
    fmpz_set_str(fmpq_numref(ans), "295147905179352825731", 10);
    fmpz_set_str(fmpq_denref(ans), "498062089990157893421", 10);
    for (i = 0; i < 1000; i++)
        fmpq_next_signed_calkin_wilf(r, r);

    if (!fmpq_equal(r, ans))
    {
        flint_printf("FAIL: signed enum from 2^64\n");
        fmpq_print(r);
        flint_printf("\n");
        fmpq_print(ans);
        flint_printf("\n");
        fflush(stdout);
        flint_abort();
    }

    fmpq_clear(r);
    fmpq_clear(ans);

    TEST_FUNCTION_END(state);
}
