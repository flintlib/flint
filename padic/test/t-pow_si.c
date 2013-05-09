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

    Copyright (C) 2011, 2012 Sebastian Pancratz

******************************************************************************/

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
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        padic_t a, b, c;
        len_t e;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);

        padic_randtest(a, state, ctx);
        padic_set(b, a, ctx);

        e = z_randint(state, 20);

        padic_pow_si(c, b, e, ctx);
        padic_pow_si(b, b, e, ctx);

        result = (padic_equal(b, c));
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("c = "), padic_print(c, ctx), printf("\n");
            printf("e = %ld\n", e);
            abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Compare with multiplication for e >= 0 and val(a) >= 0 */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        len_t N;
        padic_ctx_t ctx;

        padic_t a, b, c;
        len_t j, e;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);

        padic_randtest_int(a, state, ctx);
        padic_val(a) = FLINT_ABS(padic_val(a));

        e = n_randint(state, 50);

        padic_pow_si(b, a, e, ctx);
        padic_one(c);
        for (j = 0; j < e; j++)
            padic_mul(c, c, a, ctx);

        result = (padic_equal(b, c));
        if (!result)
        {
            printf("FAIL (cmp with multiplication):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("c = "), padic_print(c, ctx), printf("\n");
            printf("e = %ld\n", e);
            abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Compare with fmpq, check precision */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        len_t N, N2;
        padic_ctx_t ctx;

        padic_t a, b, c;
        fmpq_t s, t;
        len_t e;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);
        fmpq_init(s);
        fmpq_init(t);

        padic_randtest(a, state, ctx);
        e = n_randint(state, 50) + 1;

        N2 = N + (e - 1) * padic_val(a);

        padic_prec(b) = N2;
        padic_prec(c) = N2;
        padic_pow_si(b, a, e, ctx);

        padic_get_fmpq(s, a, ctx);
        fmpq_pow_si(t, s, e);
        padic_set_fmpq(c, t, ctx);

        result = (padic_equal(b, c));
        if (!result)
        {
            printf("FAIL (cmp with fmpq):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("c = "), padic_print(c, ctx), printf("\n");
            printf("s = "), fmpq_print(s), printf("\n");
            printf("t = "), fmpq_print(t), printf("\n");
            printf("e = %ld\n", e);
            abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);
        fmpq_clear(s);
        fmpq_clear(t);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check precision */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        len_t N_lo, N_hi, N_res;
        padic_ctx_t ctx;

        padic_t a, b, c, d, t;
        len_t e;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N_lo = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
               + PADIC_TEST_PREC_MIN;
        N_hi = N_lo + n_randint(state, 20);
        padic_ctx_init(ctx, p, FLINT_MAX(0, N_lo-10), FLINT_MAX(0, N_hi+10), PADIC_SERIES);

        padic_init2(a, N_hi);
        padic_init2(b, N_hi);
        padic_init2(t, N_lo);

        padic_randtest(a, state, ctx);
        padic_randtest(b, state, ctx);
        padic_sub(t, a, b, ctx);
        padic_add(b, b, t, ctx);
        e = n_randint(state, 50) + 1;

        N_res = N_lo + (e - 1) * FLINT_MIN(padic_val(a), padic_val(b));

        padic_init2(c, N_res);
        padic_init2(d, N_res);

        padic_pow_si(c, a, e, ctx);
        padic_pow_si(d, b, e, ctx);

        result = (padic_equal(c, d));
        if (!result)
        {
            printf("FAIL (cmp with fmpq):\n\n");
            printf("a = "), padic_print(a, ctx), printf("\n");
            printf("b = "), padic_print(b, ctx), printf("\n");
            printf("c = "), padic_print(c, ctx), printf("\n");
            printf("d = "), padic_print(d, ctx), printf("\n");
            printf("e = %ld\n", e);
            printf("N_lo  = %ld\n", N_lo);
            printf("N_hi  = %ld\n", N_hi);
            printf("N_res = %ld\n", N_res);
            abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);
        padic_clear(d);
        padic_clear(t);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

