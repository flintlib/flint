/*
      This file is part of FLINT.

      FLINT is free software: you can redistribute it and/or modify it under
      the terms of the GNU Lesser General Public License (LGPL) as published
      by the Free Software Foundation; either version 2.1 of the License, or
      (at your option) any later version.  See <https://www.gnu.org/licenses/>.
 */

/******************************************************************************

    Authored 2016 by Daniel S. Roche; US Government work in the public domain

******************************************************************************/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

/* checks that poly actually annihilates the given sequence. */
int check(const fmpz_mod_poly_t poly, const fmpz* seq, slong len,
                                                      const fmpz_mod_ctx_t ctx)
{
    fmpz_t sum, temp;
    slong d = fmpz_mod_poly_degree(poly, ctx);
    int i, j;

    if (d < 0) return 0;

    fmpz_init(sum);
    fmpz_init(temp);

    for (i=0; i < len-d; ++i)
    {
        fmpz_zero(sum);
        for (j=0; j<d; ++j)
        {
            fmpz_mod_poly_get_coeff_fmpz(temp, poly, j, ctx);
            fmpz_addmul(sum, temp, seq+(i+j));
        }
        fmpz_add(sum, sum, seq+(i+d));
        fmpz_mod(sum, sum, fmpz_mod_ctx_modulus(ctx));
        if (!fmpz_is_zero(sum))
        {
            fmpz_clear(sum);
            fmpz_clear(temp);
            return 0;
        }
    }

    fmpz_clear(sum);
    fmpz_clear(temp);
    return 1;
}

TEST_FUNCTION_START(fmpz_mod_poly_minpoly, state)
{
    int i, result;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* test random sequences */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz* seq;
        slong len;
        fmpz_mod_poly_t poly1, poly2;
        fmpz_t p;
        int j;

        fmpz_init(p);
        fmpz_randprime(p, state, 100, 0);
        fmpz_mod_ctx_set_modulus(ctx, p);

        len = n_randtest(state) % UWORD(100);
        seq = _fmpz_vec_init(len);
        for (j=0; j<len; ++j) fmpz_randtest_mod(seq+j, state, p);

        fmpz_mod_poly_init(poly1, ctx);
        fmpz_mod_poly_init(poly2, ctx);

        fmpz_mod_poly_minpoly_bm(poly1, seq, len, ctx);
        fmpz_mod_poly_minpoly_hgcd(poly2, seq, len, ctx);

        if (!check(poly1, seq, len, ctx)
            || fmpz_mod_poly_degree(poly1, ctx) > fmpz_mod_poly_degree(poly2, ctx))
        {
            flint_printf("FAIL 1:\n");
            _fmpz_vec_print(seq, len); flint_printf("\n\n");
            fmpz_mod_poly_print(poly1, ctx); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        if (!check(poly2, seq, len, ctx)
            || fmpz_mod_poly_degree(poly2, ctx) > fmpz_mod_poly_degree(poly1, ctx))
        {
            flint_printf("FAIL 2:\n");
            _fmpz_vec_print(seq, len); flint_printf("\n\n");
            fmpz_mod_poly_print(poly2, ctx); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(poly1, ctx);
        fmpz_mod_poly_clear(poly2, ctx);
        fmpz_clear(p);
        _fmpz_vec_clear(seq, len);
    }

    /* test sequences with a known generator */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz* seq;
        slong len, d;
        fmpz_mod_poly_t poly1, poly2, gen, rem;
        fmpz_t p, temp;
        int j, k;

        fmpz_init(p);
        fmpz_init(temp);
        fmpz_randprime(p, state, 100, 0);
        fmpz_mod_ctx_set_modulus(ctx, p);

        len = n_randtest(state) % UWORD(200) + 2;
        seq = _fmpz_vec_init(len);

        fmpz_mod_poly_init(poly1, ctx);
        fmpz_mod_poly_init(poly2, ctx);
        fmpz_mod_poly_init(gen, ctx);
        fmpz_mod_poly_init(rem, ctx);
        fmpz_mod_poly_randtest_monic(gen, state, n_randint(state, len/2)+2, ctx);
        d = fmpz_mod_poly_degree(gen, ctx);
        FLINT_ASSERT (d > 0);

        for (j=0; j<d; ++j) fmpz_randtest_mod(seq+j, state, p);

        for (; j<len; ++j)
        {
            fmpz_zero(seq+j);
            for (k=0; k<d; ++k)
            {
                fmpz_mod_poly_get_coeff_fmpz(temp, gen, k, ctx);
                fmpz_submul(seq+j, temp, seq + (j-d+k));
            }
            fmpz_mod(seq+j, seq+j, p);
        }
        FLINT_ASSERT(check(gen, seq, len, ctx));

        fmpz_mod_poly_minpoly_bm(poly1, seq, len, ctx);

        if (!check(poly1, seq, len, ctx))
        {
            flint_printf("FAIL 3:\n");
            _fmpz_vec_print(seq, len); flint_printf("\n\n");
            fmpz_mod_poly_print(poly1, ctx); flint_printf("\n\n");
            fmpz_mod_poly_print(gen, ctx);  flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        result = fmpz_mod_poly_degree(poly1, ctx) <= fmpz_mod_poly_degree(gen, ctx);
        if (result && fmpz_mod_poly_degree(gen, ctx) <= len/2)
        {
            fmpz_mod_poly_rem(rem, gen, poly1, ctx);
            result = fmpz_mod_poly_is_zero(rem, ctx);
        }

        if (!result)
        {
            flint_printf("FAIL 4:\n");
            _fmpz_vec_print(seq, len); flint_printf("\n\n");
            fmpz_mod_poly_print(poly1, ctx); flint_printf("\n\n");
            fmpz_mod_poly_print(gen, ctx);  flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_minpoly_hgcd(poly2, seq, len, ctx);

        if (!fmpz_mod_poly_equal(poly1, poly2, ctx))
        {
            flint_printf("FAIL 5:\n");
            _fmpz_vec_print(seq, len); flint_printf("\n\n");
            fmpz_mod_poly_print(poly2, ctx); flint_printf("\n\n");
            fmpz_mod_poly_print(gen, ctx);  flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(seq, len);
        fmpz_clear(p);
        fmpz_clear(temp);
        fmpz_mod_poly_clear(poly1, ctx);
        fmpz_mod_poly_clear(poly2, ctx);
        fmpz_mod_poly_clear(gen, ctx);
        fmpz_mod_poly_clear(rem, ctx);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
