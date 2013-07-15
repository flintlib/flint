/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_factor.h"
#include "ulong_extras.h"

void check(fmpz_t n)
{
    fmpz_factor_t factor;
    fmpz_t m;

    fmpz_factor_init(factor);
    fmpz_init(m);

    fmpz_factor(factor, n);
    fmpz_factor_expand(m, factor);

    if (!fmpz_equal(n, m))
    {
        printf("ERROR: factors do not unfactor to original number!\n");

        printf("input: ");
        fmpz_print(n);
        printf("\n");

        printf("computed factors: ");
        fmpz_factor_print(factor);
        printf("\n");

        printf("value: ");
        fmpz_print(m);
        printf("\n");

        abort();
    }

    fmpz_clear(m);
    fmpz_factor_clear(factor);
}

int main(void)
{
    int i, j;
    fmpz_t x;
    mpz_t y;

    printf("factor....");
    fflush(stdout);

    fmpz_init(x);
    mpz_init(y);

    /* Some corner cases */
    fmpz_set_ui(x, ULONG_MAX);
    check(x);
    fmpz_set_si(x, LONG_MAX);
    check(x);
    fmpz_set_si(x, LONG_MIN);
    check(x);
    fmpz_set_si(x, COEFF_MAX);
    check(x);
    fmpz_set_si(x, COEFF_MIN);
    check(x);

    /* Small integers */
    for (i = -10000; i < 10000; i++)
    {
        fmpz_set_si(x, i);
        check(x);
    }

    /* Powers */
    for (i = 1; i < 250; i++)
    {
        for (j = 0; j < 250; j++)
        {
            fmpz_set_ui(x, i);
            fmpz_pow_ui(x, x, j);
            check(x);
        }
    }

    /* Factorials */
    for (i = 0; i < 1000; i++)
    {
        mpz_fac_ui(y, i);
        fmpz_set_mpz(x, y);
        check(x);
    }

    /* Powers of factorials */
    for (i = 0; i < 100; i++)
    {
        for (j = 1; j < 5; j++)
        {
            mpz_fac_ui(y, i);
            fmpz_set_mpz(x, y);
            fmpz_pow_ui(x, x, j);
            check(x);
        }
    }

    /* Large negative integers */
    fmpz_set_ui(x, 10);
    fmpz_pow_ui(x, x, 100);
    fmpz_neg(x, x);
    check(x);
    mpz_fac_ui(y, 50);
    mpz_neg(y, y);
    fmpz_set_mpz(x, y);
    check(x);

    fmpz_clear(x);
    mpz_clear(y);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
