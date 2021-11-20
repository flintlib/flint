/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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
    FLINT_TEST_INIT(state);

    flint_printf("fdiv_qr_preinvn....");
    fflush(stdout);

    

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, r;
        mpz_t d, e, f, g, h, s;
        fmpz_preinvn_t inv;

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

        fmpz_randtest(a, state, 400);
        fmpz_randtest_not_zero(b, state, 400);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);

        fmpz_preinvn_init(inv, b);
        fmpz_fdiv_qr_preinvn(c, r, a, b, inv);
        mpz_fdiv_qr(f, s, d, e);

        fmpz_get_mpz(g, c);
        fmpz_get_mpz(h, r);

        result = (mpz_cmp(f, g) == 0 && mpz_cmp(h, s) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf
                ("d = %Zd, e = %Zd, f = %Zd, g = %Zd, h = %Zd, s = %Zd\n", d,
                 e, f, g, h, s);
            fflush(stdout);
            flint_abort();
        }

        fmpz_preinvn_clear(inv);
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
        fmpz_preinvn_t inv;

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

        fmpz_preinvn_init(inv, b);
        fmpz_fdiv_qr_preinvn(a, r, a, b, inv);
        mpz_fdiv_qr(f, s, d, e);

        fmpz_get_mpz(g, a);
        fmpz_get_mpz(h, r);

        result = (mpz_cmp(f, g) == 0 && mpz_cmp(h, s) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf
                ("d = %Zd, e = %Zd, f = %Zd, g = %Zd, h = %Zd, s = %Zd\n", d,
                 e, f, g, h, s);
            fflush(stdout);
            flint_abort();
        }

        fmpz_preinvn_clear(inv);
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
        fmpz_preinvn_t inv;

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

        fmpz_preinvn_init(inv, b);
        fmpz_fdiv_qr_preinvn(b, r, a, b, inv);
        mpz_fdiv_qr(f, s, d, e);

        fmpz_get_mpz(g, b);
        fmpz_get_mpz(h, r);

        result = (mpz_cmp(f, g) == 0 && mpz_cmp(h, s) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf
                ("d = %Zd, e = %Zd, f = %Zd, g = %Zd, h = %Zd, s = %Zd\n", d,
                 e, f, g, h, s);
            fflush(stdout);
            flint_abort();
        }

        fmpz_preinvn_clear(inv);
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
        fmpz_preinvn_t inv;

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

        fmpz_preinvn_init(inv, b);
        fmpz_fdiv_qr_preinvn(c, a, a, b, inv);
        mpz_fdiv_qr(f, s, d, e);

        fmpz_get_mpz(g, c);
        fmpz_get_mpz(h, a);

        result = (mpz_cmp(f, g) == 0 && mpz_cmp(h, s) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf
                ("d = %Zd, e = %Zd, f = %Zd, g = %Zd, h = %Zd, s = %Zd\n", d,
                 e, f, g, h, s);
            fflush(stdout);
            flint_abort();
        }

        fmpz_preinvn_clear(inv);
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
        fmpz_preinvn_t inv;

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

        fmpz_preinvn_init(inv, b);
        fmpz_fdiv_qr_preinvn(c, b, a, b, inv);
        mpz_fdiv_qr(f, s, d, e);

        fmpz_get_mpz(g, c);
        fmpz_get_mpz(h, b);

        result = (mpz_cmp(f, g) == 0 && mpz_cmp(h, s) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf
                ("d = %Zd, e = %Zd, f = %Zd, g = %Zd, h = %Zd, s = %Zd\n", d,
                 e, f, g, h, s);
            fflush(stdout);
            flint_abort();
        }

        fmpz_preinvn_clear(inv);
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

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
