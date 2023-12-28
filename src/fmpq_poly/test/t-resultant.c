/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpq_poly.h"

#ifdef __GNUC__
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Woverlength-strings"
#endif

TEST_FUNCTION_START(fmpq_poly_resultant, state)
{
    int i, result;

    /* Check res(f, g) == (-1)^(deg f deg g) res(g, f) */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g;
        fmpq_t x, y;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_init(x);
        fmpq_init(y);

        fmpq_poly_randtest(f, state, n_randint(state, 60), 60);
        fmpq_poly_randtest(g, state, n_randint(state, 60), 60);

        fmpq_poly_resultant(x, f, g);
        fmpq_poly_resultant(y, g, f);
        if ((fmpq_poly_degree(f) * fmpq_poly_degree(g)) % 2)
            fmpq_neg(y, y);

        result = fmpq_equal(x, y);
        if (!result)
        {
            flint_printf("FAIL (res(f,g) == (-1)^(m * n) res(g, f)):\n");
            flint_printf("f = "), fmpq_poly_print(f), flint_printf("\n\n");
            flint_printf("g = "), fmpq_poly_print(g), flint_printf("\n\n");
            flint_printf("x = "), fmpq_print(x), flint_printf("\n\n");
            flint_printf("y = "), fmpq_print(y), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_clear(x);
        fmpq_clear(y);
    }

    /* Check res(f h, g) == res(f, g) res(h, g) */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h;
        fmpq_t x, y, z;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        fmpq_poly_randtest(f, state, n_randint(state, 60), 60);
        fmpq_poly_randtest(g, state, n_randint(state, 60), 60);
        fmpq_poly_randtest(h, state, n_randint(state, 60), 60);

        fmpq_poly_resultant(y, f, g);
        fmpq_poly_resultant(z, h, g);
        fmpq_mul(y, y, z);
        fmpq_poly_mul(f, f, h);
        fmpq_poly_resultant(x, f, g);

        result = fmpq_equal(x, y);
        if (!result)
        {
            flint_printf("FAIL (res(f h, g) == res(f, g) res(h, g)):\n");
            flint_printf("f = "), fmpq_poly_print(f), flint_printf("\n\n");
            flint_printf("g = "), fmpq_poly_print(g), flint_printf("\n\n");
            flint_printf("h = "), fmpq_poly_print(h), flint_printf("\n\n");
            flint_printf("x = "), fmpq_print(x), flint_printf("\n\n");
            flint_printf("y = "), fmpq_print(y), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* fredrik's test case */
    {
        fmpq_poly_t f, g;
        fmpq_t x, y;
        int result;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_init(x);
        fmpq_init(y);

        fmpq_poly_set_str(f, "49  16702090503 -23810415210 7561766512 801950253"
             " 56796743 40735271 -15934 820601 -2712604160 -1577466 0 0 -7967 0"
             " 0 0 -14491973 0 6566138489 -55769 0 130523361 4071137 15934"
             " -501921 -59067338 63860755253 23924901 -15934 -262911 -7967"
             " -4389817 0 185876611072 58470888545 130523361 -63736 -130618965"
             " -39835 0 7967 0 55769 -7967 103571 111298990 47802 -3808226"
             " -3800259");

        fmpq_poly_set_str(g, "59  -458395/219902324736 151585/4581298432"
           " 112595/219902324736 -2016245/54975581184 0 35/73300774912 0"
           " -234880919/219902324736 7/219902324736 -7/1278501888"
           " -6055/109951162368 7/27487790592 -504623/73300774912"
           " 53673977/219902324736 0 611667/73300774912 -497/13743895296"
           " 0 -6265/219902324736 2446675/73300774912 2345/219902324736"
           " -371/73300774912 -427/6871947648 -3758096377/219902324736"
           " 20595995/109951162368 -256459/73300774912 0 33690223/73300774912"
           " -229369/219902324736 93205/219902324736 -7/107374182"
           " -133/219902324736 -665/13743895296 -146503/219902324736 0"
           " 7/219902324736 66633/73300774912 -855190385/219902324736"
           " 229355/219902324736 0 161/219902324736 887299/219902324736"
           " -427/7582838784 -611667/18325193728 -7/5114007552 833/54975581184"
           " -7/109951162368 -5402264413/219902324736 7/5114007552 35/9162596864"
           " 1133545/219902324736 -151319/73300774912 0 7/219902324736"
           " 7/54975581184 0 -10367/109951162368 7/54975581184 -161/109951162368");

        fmpq_poly_resultant(x, f, g);
        fmpq_poly_resultant(y, g, f);

        if ((fmpq_poly_degree(f) * fmpq_poly_degree(g)) % 2)
            fmpq_neg(y, y);

        result = fmpq_equal(x, y);
        if (!result)
        {
            flint_printf("FAIL (res(f,g) == (-1)^(m * n) res(g, f)):\n");
            flint_printf("x = "), fmpq_print(x), flint_printf("\n\n");
            flint_printf("y = "), fmpq_print(y), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_clear(x);
        fmpq_clear(y);
    }

    TEST_FUNCTION_END(state);
}
#ifdef __GNUC__
# pragma GCC diagnostic pop
#endif
