/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
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
    fmpq_set_si(ans, 191, 1287);
    for (i = 0; i < 1000000; i++)
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
    fmpq_set_si(ans, -191, 1096);
    for (i = 0; i < 1000000; i++)
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
    fmpz_set_str(fmpq_numref(ans), "18538977794078099355206", 10);
    fmpz_set_str(fmpq_denref(ans), "22062305912156623710275", 10);
    for (i = 0; i < 1000000; i++)
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
    fmpz_set_str(fmpq_numref(ans), "15015649675999575000951", 10);
    fmpz_set_str(fmpq_denref(ans), "18538977794078099356211", 10);
    for (i = 0; i < 1000000; i++)
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
