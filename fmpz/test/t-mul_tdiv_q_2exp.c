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

    Copyright (C) 2009 William Hart
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("mul_tdiv_q_2exp....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c;
        mpz_t d, e, f, g;
        ulong exp;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest(a, state, 200);
        fmpz_randtest(b, state, 200);
        exp = n_randint(state, 200);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        fmpz_mul_tdiv_q_2exp(c, a, b, exp);
        mpz_mul(f, d, e);
        mpz_tdiv_q_2exp(f, f, exp);

        fmpz_get_mpz(g, c);

        result = (mpz_cmp(f, g) == 0);
        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, g = %Zd, exp = %lu\n",
                d, e, f, g, exp);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
    }

    /* Check aliasing of a and b */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, c;
        mpz_t d, f, g;
        ulong exp;

        fmpz_init(a);
        fmpz_init(c);

        mpz_init(d);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest(a, state, 200);
        exp = n_randint(state, 200);
        fmpz_get_mpz(d, a);

        fmpz_mul_tdiv_q_2exp(c, a, a, exp);
        mpz_mul(f, d, d);
        mpz_tdiv_q_2exp(f, f, exp);

        fmpz_get_mpz(g, c);

        result = (mpz_cmp(f, g) == 0);
        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf("d = %Zd, f = %Zd, g = %Zd, exp = %lu\n", d, f, g, exp);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);

        mpz_clear(d);
        mpz_clear(f);
        mpz_clear(g);
    }

    /* Test aliasing of a and c */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t d, e, f, g;
        ulong exp;

        fmpz_init(a);
        fmpz_init(b);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest(a, state, 200);
        fmpz_randtest(b, state, 200);
        exp = n_randint(state, 200);
        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        fmpz_mul_tdiv_q_2exp(a, a, b, exp);
        mpz_mul(f, d, e);
        mpz_tdiv_q_2exp(f, f, exp);

        fmpz_get_mpz(g, a);

        result = (mpz_cmp(f, g) == 0);
        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, g = %Zd, exp = %lu\n",
                d, e, f, g, exp);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
    }

    /* Test aliasing of b and c */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t d, e, f, g;
        ulong exp;

        fmpz_init(a);
        fmpz_init(b);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);

        fmpz_randtest(a, state, 200);
        fmpz_randtest(b, state, 200);
        exp = n_randint(state, 200);
        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        fmpz_mul_tdiv_q_2exp(b, a, b, exp);
        mpz_mul(f, d, e);
        mpz_tdiv_q_2exp(f, f, exp);

        fmpz_get_mpz(g, b);

        result = (mpz_cmp(f, g) == 0);
        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, g = %Zd, exp = %lu\n",
                d, e, f, g, exp);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
