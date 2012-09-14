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
    Copyright (C) 2012 Andres Goens

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "fq.h"
#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("frobenius... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long d;
        fq_ctx_t ctx;

        fq_t a, b, c;
        long e;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");

        fq_init(a);
        fq_init(b);
        fq_init(c);

        fq_randtest(a, state, ctx);
        fq_set(b, a);
        e = n_randint(state, 10) % d;

        fq_frobenius(c, b, e, ctx);
        fq_frobenius(b, b, e, ctx);

        result = (fq_equal(b, c));
        if (!result)
        {
            printf("FAIL (alias):\n\n");
            printf("a = "), fq_print_pretty(a, ctx), printf("\n");
            printf("b = "), fq_print_pretty(b, ctx), printf("\n");
            printf("c = "), fq_print_pretty(c, ctx), printf("\n");
            printf("e = %ld\n", e);
            abort();
        }

        fq_clear(a);
        fq_clear(b);
        fq_clear(c);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }

    /* Check sigma^e(x) == x^{p^e}  */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long d;
        fq_ctx_t ctx;

        fq_t a, b, c;
        long e;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");

        fq_init(a);
        fq_init(b);
        fq_init(c);

        fq_randtest(a, state, ctx);
        e = n_randint(state, 10) % d;

        fq_frobenius(b, a, e, ctx);
        {
            fmpz_t t;

            fmpz_init(t);
            fmpz_pow_ui(t, p, e);
            fq_pow(c, a, t, ctx);
            fmpz_clear(t);
        }

        result = (fq_equal(b,c));
        if (!result)
        {
            printf("FAIL (sigma^e(x) = x^{p^e}):\n\n");
            printf("a = "), fq_print_pretty(a, ctx), printf("\n");
            printf("b = "), fq_print_pretty(b, ctx), printf("\n");
            printf("c = "), fq_print_pretty(c, ctx), printf("\n");
            printf("e = %ld\n", e);
            abort();
        }

        fq_clear(a);
        fq_clear(b);
        fq_clear(c);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }

    /* Check sigma^e(x + y) = sigma^e(x) + sigma^e(y) */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long d;
        fq_ctx_t ctx;

        fq_t a, b, s, s1, s2, lhs, rhs;
        long e;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");

        fq_init(a);
        fq_init(b);
        fq_init(s);
        fq_init(s1);
        fq_init(s2);
        fq_init(lhs);
        fq_init(rhs);

        fq_randtest(a, state, ctx);
        fq_randtest(b, state, ctx);
        e = n_randint(state, 10) % d;

        fq_add(s, a, b, ctx);
        fq_frobenius(lhs, s, e, ctx);
        fq_frobenius(s1, a, e, ctx);
        fq_frobenius(s2, b, e, ctx);
        fq_add(rhs, s1, s2, ctx);

        result = (fq_equal(lhs, rhs));
        if (!result)
        {
            printf("FAIL (sigma(a+b) = sigma(a) + sigma(b)):\n\n");
            printf("a = "), fq_print_pretty(a, ctx), printf("\n");
            printf("b = "), fq_print_pretty(b, ctx), printf("\n");
            printf("s = "), fq_print_pretty(s, ctx), printf("\n");
            printf("s1 = "), fq_print_pretty(s1, ctx), printf("\n");
            printf("s2 = "), fq_print_pretty(s2, ctx), printf("\n");
            printf("lhs = "), fq_print_pretty(lhs, ctx), printf("\n");
            printf("rhs = "), fq_print_pretty(rhs, ctx), printf("\n");
            printf("e = %ld\n", e);
            abort();
        }

        fq_clear(a);
        fq_clear(b);
        fq_clear(s);
        fq_clear(s1);
        fq_clear(s2);
        fq_clear(lhs);
        fq_clear(rhs);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }

    /* Check sigma^e(x * y) = sigma^e(x) * sigma^e(y) on Zq */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long d;
        fq_ctx_t ctx;

        fq_t a, b, s, s1, s2, lhs, rhs;
        long e;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        fq_ctx_init_conway(ctx, p, d, "a");

        fq_init(a);
        fq_init(b);
        fq_init(s);
        fq_init(s1);
        fq_init(s2);
        fq_init(lhs);
        fq_init(rhs);

        fq_randtest(a, state, ctx);
        fq_randtest(b, state, ctx);
        e = n_randint(state, 10) % d;

        fq_mul(s, a, b, ctx);
        fq_frobenius(lhs, s, e, ctx);
        fq_frobenius(s1, a, e, ctx);
        fq_frobenius(s2, b, e, ctx);
        fq_mul(rhs, s1, s2, ctx);

        result = (fq_equal(lhs, rhs));
        if (!result)
        {
            printf("FAIL (sigma(a*b) = sigma(a) * sigma(b)):\n\n");
            printf("a = "), fq_print_pretty(a, ctx), printf("\n");
            printf("b = "), fq_print_pretty(b, ctx), printf("\n");
            printf("s = "), fq_print_pretty(s, ctx), printf("\n");
            printf("s1 = "), fq_print_pretty(s1, ctx), printf("\n");
            printf("s2 = "), fq_print_pretty(s2, ctx), printf("\n");
            printf("lhs = "), fq_print_pretty(lhs, ctx), printf("\n");
            printf("rhs = "), fq_print_pretty(rhs, ctx), printf("\n");
            printf("e = %ld\n", e);
            abort();
        }

        fq_clear(a);
        fq_clear(b);
        fq_clear(s);
        fq_clear(s1);
        fq_clear(s2);
        fq_clear(lhs);
        fq_clear(rhs);

        fmpz_clear(p);
        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

