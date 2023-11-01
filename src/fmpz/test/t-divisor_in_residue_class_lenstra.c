/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_divisor_in_residue_class_lenstra, state)
{
    int i, result;

    /* test factors of composites are found */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p, a, r, s, d, g;
        int res;

        fmpz_init(p);
        fmpz_init(a);
        fmpz_init(d);
        fmpz_init(r);
        fmpz_init(s);
        fmpz_init(g);

        do {
           do {
              fmpz_randbits(p, state, n_randint(state, 100) + 2);
           } while (fmpz_cmp_ui(p, 2) < 0);
           do {
              fmpz_randbits(a, state, n_randint(state, 100) + 2);
           } while (fmpz_cmp_ui(a, 2) < 0);

           fmpz_mul(p, p, a);

           fmpz_root(s, p, 3); /* cube root of p */
           fmpz_randbits(r, state, (2*fmpz_bits(p))/3);
           fmpz_abs(r, r);

           fmpz_mul(s, s, r); /* s now between cube root and p */

           fmpz_mod(r, a, s);

           fmpz_gcd(g, r, s);
        } while (!fmpz_is_one(g));

        result = ((res = fmpz_divisor_in_residue_class_lenstra(d, p, r, s)) && !fmpz_is_one(d)
           && !fmpz_equal(d, p) && fmpz_divisible(p, d));
        if (!result)
        {
            flint_printf("FAIL:\n");
            printf("%d\n", res);
            fmpz_print(p); printf("\n");
            fmpz_print(r); printf("\n");
            fmpz_print(s); printf("\n");
            fmpz_print(a); printf("\n");
            fmpz_print(d); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_clear(a);
        fmpz_clear(r);
        fmpz_clear(s);
        fmpz_clear(g);
        fmpz_clear(d);
    }

    /* test factors of primes are not found */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p, r, s, d, g;
        int res;

        fmpz_init(p);
        fmpz_init(d);
        fmpz_init(r);
        fmpz_init(s);
        fmpz_init(g);

        do {
           do {
              fmpz_randbits(p, state, n_randint(state, 100) + 2);
           } while (!fmpz_is_probabprime_BPSW(p));

           fmpz_root(s, p, 3); /* cube root of p */
           fmpz_randbits(r, state, (2*fmpz_bits(p))/3);
           fmpz_abs(r, r);

           fmpz_mul(s, s, r); /* s now between cube root and p */

           fmpz_randm(r, state, s);

           fmpz_gcd(g, r, s);
        } while (!fmpz_is_one(g) || fmpz_is_one(s));

        result = (!(res = fmpz_divisor_in_residue_class_lenstra(d, p, r, s)));
        if (!result)
        {
            flint_printf("FAIL:\n");
            printf("%d\n", res);
            fmpz_print(p); printf("\n");
            fmpz_print(r); printf("\n");
            fmpz_print(s); printf("\n");
            fmpz_print(d); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_clear(r);
        fmpz_clear(s);
        fmpz_clear(g);
        fmpz_clear(d);
    }

    TEST_FUNCTION_END(state);
}
