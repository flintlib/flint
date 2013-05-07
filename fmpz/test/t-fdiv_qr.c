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

    printf("fdiv_qr....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, r;
        mpz_t d, e, f, g, h, s;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(r);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);
        mpz_init(h);
        mpz_init(s);

        fmpz_randtest(a, state, 200);
        fmpz_randtest_not_zero(b, state, 200);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        fmpz_fdiv_qr(c, r, a, b);
        mpz_fdiv_qr(f, s, d, e);

        fmpz_get_mpz(g, c);
        fmpz_get_mpz(h, r);

        result = (mpz_cmp(f, g) == 0 && mpz_cmp(h, s) == 0);
        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf
                ("d = %Zd, e = %Zd, f = %Zd, g = %Zd, h = %Zd, s = %Zd\n", d,
                 e, f, g, h, s);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(r);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
        mpz_clear(h);
        mpz_clear(s);
    }

    /* Check aliasing of c and a */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, r;
        mpz_t d, e, f, g, h, s;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(r);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);
        mpz_init(h);
        mpz_init(s);

        fmpz_randtest(a, state, 200);
        fmpz_randtest_not_zero(b, state, 200);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        fmpz_fdiv_qr(a, r, a, b);
        mpz_fdiv_qr(f, s, d, e);

        fmpz_get_mpz(g, a);
        fmpz_get_mpz(h, r);

        result = (mpz_cmp(f, g) == 0 && mpz_cmp(h, s) == 0);
        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf
                ("d = %Zd, e = %Zd, f = %Zd, g = %Zd, h = %Zd, s = %Zd\n", d,
                 e, f, g, h, s);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(r);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
        mpz_clear(h);
        mpz_clear(s);
    }

    /* Check aliasing of c and b */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, r;
        mpz_t d, e, f, g, h, s;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(r);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);
        mpz_init(h);
        mpz_init(s);

        fmpz_randtest(a, state, 200);
        fmpz_randtest_not_zero(b, state, 200);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        fmpz_fdiv_qr(b, r, a, b);
        mpz_fdiv_qr(f, s, d, e);

        fmpz_get_mpz(g, b);
        fmpz_get_mpz(h, r);

        result = (mpz_cmp(f, g) == 0 && mpz_cmp(h, s) == 0);
        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf
                ("d = %Zd, e = %Zd, f = %Zd, g = %Zd, h = %Zd, s = %Zd\n", d,
                 e, f, g, h, s);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(r);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
        mpz_clear(h);
        mpz_clear(s);
    }

    /* Check aliasing of r and a */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, r;
        mpz_t d, e, f, g, h, s;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(r);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);
        mpz_init(h);
        mpz_init(s);

        fmpz_randtest(a, state, 200);
        fmpz_randtest_not_zero(b, state, 200);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        fmpz_fdiv_qr(c, a, a, b);
        mpz_fdiv_qr(f, s, d, e);

        fmpz_get_mpz(g, c);
        fmpz_get_mpz(h, a);

        result = (mpz_cmp(f, g) == 0 && mpz_cmp(h, s) == 0);
        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf
                ("d = %Zd, e = %Zd, f = %Zd, g = %Zd, h = %Zd, s = %Zd\n", d,
                 e, f, g, h, s);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(r);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
        mpz_clear(h);
        mpz_clear(s);
    }

    /* Check aliasing of r and b */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, r;
        mpz_t d, e, f, g, h, s;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(r);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(g);
        mpz_init(h);
        mpz_init(s);

        fmpz_randtest(a, state, 200);
        fmpz_randtest_not_zero(b, state, 200);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        fmpz_fdiv_qr(c, b, a, b);
        mpz_fdiv_qr(f, s, d, e);

        fmpz_get_mpz(g, c);
        fmpz_get_mpz(h, b);

        result = (mpz_cmp(f, g) == 0 && mpz_cmp(h, s) == 0);
        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf
                ("d = %Zd, e = %Zd, f = %Zd, g = %Zd, h = %Zd, s = %Zd\n", d,
                 e, f, g, h, s);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(r);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(g);
        mpz_clear(h);
        mpz_clear(s);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
