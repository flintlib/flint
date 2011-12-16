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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpq.h"
#include "padic.h"
#include "long_extras.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("pow_si... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        padic_t a, b, c;
        long e;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 100) + 1;
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(b, ctx);
        padic_init(c, ctx);

        padic_randtest(a, state, ctx);
        padic_set(b, a, ctx);

        e = z_randint(state, 20);

        padic_pow_si(c, b, e, ctx);
        padic_pow_si(b, b, e, ctx);

        result = (padic_equal(b, c, ctx));
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("c = "), padic_print(c, ctx), printf("\n");
            printf("e = %ld\n", e);
            abort();
        }

        padic_clear(a, ctx);
        padic_clear(b, ctx);
        padic_clear(c, ctx);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Compare with multiplication for e > 0 and val(a) >= 0 */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx;

        padic_t a, b, c;
        long j, e;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 100);
        padic_ctx_init(ctx, p, N, PADIC_SERIES);

        padic_init(a, ctx);
        padic_init(b, ctx);
        padic_init(c, ctx);

        padic_randtest(a, state, ctx);
        padic_val(a) = FLINT_ABS(padic_val(a));

        e = n_randint(state, 50) + 1;

        padic_pow_si(b, a, e, ctx);
        _padic_one(c);
        for (j = 0; j < e; j++)
            padic_mul(c, c, a, ctx);

        result = (padic_equal(b, c, ctx));
        if (!result)
        {
            printf("FAIL (cmp with multiplication):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("c = "), padic_print(c, ctx), printf("\n");
            printf("e = %ld\n", e);
            abort();
        }

        padic_clear(a, ctx);
        padic_clear(b, ctx);
        padic_clear(c, ctx);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Compare with fmpq, check precision */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N;
        padic_ctx_t ctx1, ctx2;

        padic_t a, b, c;
        fmpq_t s, t;
        long e;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N = z_randint(state, 100);
        padic_ctx_init(ctx1, p, N, PADIC_SERIES);

        _padic_init(a);
        _padic_init(b);
        _padic_init(c);
        fmpq_init(s);
        fmpq_init(t);

        padic_randtest(a, state, ctx1);
        e = n_randint(state, 50) + 1;

        padic_ctx_init(ctx2, p, N + (e - 1) * padic_val(a), PADIC_SERIES);

        padic_pow_si(b, a, e, ctx2);

        padic_get_fmpq(s, a, ctx1);
        fmpq_pow_si(t, s, e);
        padic_set_fmpq(c, t, ctx2);

        result = (padic_equal(b, c, ctx2));
        if (!result)
        {
            printf("FAIL (cmp with fmpq):\n\n");
            printf("a = "), padic_print(a, ctx1), printf("\n");
            printf("b = "), padic_print(b, ctx2), printf("\n");
            printf("c = "), padic_print(c, ctx2), printf("\n");
            printf("s = "), fmpq_print(s), printf("\n");
            printf("t = "), fmpq_print(t), printf("\n");
            printf("e = %ld\n", e);
            abort();
        }

        _padic_clear(a);
        _padic_clear(b);
        _padic_clear(c);
        fmpq_clear(s);
        fmpq_clear(t);

        fmpz_clear(p);
        padic_ctx_clear(ctx1);
        padic_ctx_clear(ctx2);
    }

    /* Check precision */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t p;
        long N_lo, N_hi, N_res;
        padic_ctx_t ctx_lo, ctx_hi, ctx_res;

        padic_t a, b, c, d, t;
        long e;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 5, 1));
        N_lo = z_randint(state, 100);
        N_hi = N_lo + n_randint(state, 20);
        padic_ctx_init(ctx_lo, p, N_lo, PADIC_SERIES);
        padic_ctx_init(ctx_hi, p, N_hi, PADIC_SERIES);

        _padic_init(a);
        _padic_init(b);
        _padic_init(c);
        _padic_init(d);
        _padic_init(t);

        padic_randtest(a, state, ctx_hi);
        padic_randtest(b, state, ctx_hi);
        padic_sub(t, a, b, ctx_hi);
        padic_reduce(t, ctx_lo);
        padic_add(b, b, t, ctx_hi);
        e = n_randint(state, 50) + 1;

        N_res = N_lo + (e - 1) * FLINT_MIN(padic_val(a), padic_val(b));
        padic_ctx_init(ctx_res, p, N_res, PADIC_SERIES);

        padic_pow_si(c, a, e, ctx_res);
        padic_pow_si(d, b, e, ctx_res);

        result = (padic_equal(c, d, ctx_res));
        if (!result)
        {
            printf("FAIL (cmp with fmpq):\n\n");
            printf("a = "), padic_print(a, ctx_hi), printf("\n");
            printf("b = "), padic_print(b, ctx_hi), printf("\n");
            printf("c = "), padic_print(c, ctx_res), printf("\n");
            printf("d = "), padic_print(d, ctx_res), printf("\n");
            printf("e = %ld\n", e);
            printf("N_lo  = %ld\n", N_lo);
            printf("N_hi  = %ld\n", N_hi);
            printf("N_res = %ld\n", N_res);
            abort();
        }

        _padic_clear(a);
        _padic_clear(b);
        _padic_clear(c);
        _padic_clear(d);
        _padic_clear(t);

        fmpz_clear(p);
        padic_ctx_clear(ctx_lo);
        padic_ctx_clear(ctx_hi);
        padic_ctx_clear(ctx_res);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

