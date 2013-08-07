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

    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include "ulong_extras.h"
#include "padic.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("padic_val_fac... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing for padic_val_fac() */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t a, b, c, p;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(p);

        fmpz_randtest_unsigned(a, state, (mp_bitcnt_t) (1.5 * FLINT_BITS));
        fmpz_set(b, a);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        padic_val_fac(c, b, p);
        padic_val_fac(b, b, p);

        result = fmpz_equal(b, c);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("b = "), fmpz_print(b), printf("\n");
            printf("c = "), fmpz_print(c), printf("\n");
            printf("p = "), fmpz_print(p), printf("\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(p);
    }

    /* Check correctness for padic_val_fac_ui(), p == 2 */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t a, b, p;

        ulong s, t, N;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(p);

        N = n_randint(state, 1L < 13);
        fmpz_set_ui(p, 2);
        fmpz_fac_ui(a, N);

        s = padic_val_fac_ui_2(N);
        t = fmpz_remove(b, a, p);

        result = (s == t);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("b = "), fmpz_print(b), printf("\n");
            printf("N = %lu\n", N);
            printf("s = %lu\n", s);
            printf("t = %lu\n", t);
            printf("p = "), fmpz_print(p), printf("\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(p);
    }

    /* Check correctness for padic_val_fac_ui(), any p */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t a, b, p;

        ulong s, t, N;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(p);

        N = n_randint(state, 1L < 13);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, FLINT_BITS - 4), 0));
        fmpz_fac_ui(a, N);

        s = padic_val_fac_ui(N, p);
        t = fmpz_remove(b, a, p);

        result = (s == t);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("b = "), fmpz_print(b), printf("\n");
            printf("N = %lu\n", N);
            printf("s = %lu\n", s);
            printf("t = %lu\n", t);
            printf("p = "), fmpz_print(p), printf("\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(p);
    }

    /* Compare padic_val_fac_ui() with padic_val_fac() */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t a, p, t, z;

        ulong s, n;

        fmpz_init(a);
        fmpz_init(p);
        fmpz_init(t);
        fmpz_init(z);

        n = n_randint(state, 1L < 13);
        fmpz_set_ui(z, n);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, FLINT_BITS - 4), 0));
        fmpz_fac_ui(a, n);

        s = padic_val_fac_ui(n, p);
        padic_val_fac(t, z, p);

        result = (fmpz_equal_ui(t, s));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), fmpz_print(a), printf("\n");
            printf("n = %lu\n", n);
            printf("s = %lu\n", s);
            printf("t = "), fmpz_print(t), printf("\n");
            printf("p = "), fmpz_print(p), printf("\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(p);
        fmpz_clear(t);
        fmpz_clear(z);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}
