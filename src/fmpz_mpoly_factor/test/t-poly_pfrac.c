/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"
#include "fmpz_mpoly_factor.h"
#include "fmpq_poly.h"

void _test_pfrac(
    fmpz_poly_struct * c,
    fmpz_poly_pfrac_t I,
    const fmpz_poly_struct * b,
    slong r,
    flint_rand_t state)
{
    slong i, j;
    int success, found_bad;
    fmpq_poly_struct * cQ, * bQ, * prod_bQ, * inv_prod_bQ;
    fmpq_poly_t aQ, pQ, G, S;
    fmpz_poly_t a, t, t1;

    cQ = FLINT_ARRAY_ALLOC(r, fmpq_poly_struct);
    bQ = FLINT_ARRAY_ALLOC(r, fmpq_poly_struct);
    prod_bQ = FLINT_ARRAY_ALLOC(r, fmpq_poly_struct);
    inv_prod_bQ = FLINT_ARRAY_ALLOC(r, fmpq_poly_struct);
    for (i = 0; i < r; i++)
    {
        fmpq_poly_init(cQ + i);
        fmpq_poly_init(bQ + i);
        fmpq_poly_init(prod_bQ + i);
        fmpq_poly_init(inv_prod_bQ + i);
    }

    fmpz_poly_init(a);
    fmpz_poly_init(t);
    fmpz_poly_init(t1);
    fmpq_poly_init(aQ);
    fmpq_poly_init(pQ);
    fmpq_poly_init(G);
    fmpq_poly_init(S);

    fmpq_poly_one(pQ);

    for (i = 0; i < r; i++)
    {
        fmpq_poly_set_fmpz_poly(bQ + i, b + i);
        fmpq_poly_mul(pQ, pQ, bQ + i);
    }

    success = 1;
    for (i = 0; i < r; i++)
    {
        fmpq_poly_divrem(prod_bQ + i, G, pQ, bQ + i);
        fmpq_poly_xgcd(G, S, inv_prod_bQ + i, bQ + i, prod_bQ + i);
        if (!fmpq_poly_is_one(G))
            success = 0;
    }

    if (success != fmpz_poly_pfrac_precompute(I, b, r))
    {
        flint_printf("FAIL: check precompute\n");
        fflush(stdout);
        flint_abort();
    }

    if (!success)
        goto cleanup;

    for (j = 0; j < 20; j++)
    {
        if (j % 6)
        {
            fmpz_poly_zero(a);
            fmpz_poly_one(t);
            for (i = 0; i < r; i++)
            {
                fmpz_poly_randtest(t1, state, fmpz_poly_degree(b + i),
                                                    2 + n_randint(state, 200));
                fmpz_poly_mul(t1, t1, t);
                fmpz_poly_mul(a, a, b + i);
                fmpz_poly_add(a, a, t1);
                fmpz_poly_mul(t, t, b + i);
            }
        }
        else
        {
            fmpz_poly_randtest(a, state, n_randint(state, pQ->length),
                                                    2 + n_randint(state, 200));
        }

        FLINT_ASSERT(a->length < pQ->length);

        fmpq_poly_set_fmpz_poly(aQ, a);

        found_bad = 0;
        for (i = 0; i < r; i++)
        {
            fmpq_poly_mul(S, aQ, inv_prod_bQ + i);
            fmpq_poly_rem(cQ + i, S, bQ + i);
            if (!fmpz_is_one(cQ[i].den))
                found_bad = 1;
        }

        if (fmpz_poly_pfrac_precomp(c, a, I))
        {
            if (found_bad)
            {
                flint_printf("FAIL: precomp should have failed");
                fflush(stdout);
                flint_abort();
            }

            for (i = 0; i < r; i++)
            {
                fmpq_poly_set_fmpz_poly(S, c + i);
                if (!fmpq_poly_equal(S, cQ + i))
                {
                    flint_printf("FAIL: precomp produced wrong answer\n");
                    fflush(stdout);
                    flint_abort();
                }
            }
        }
        else
        {
            if (!found_bad)
            {
                flint_printf("FAIL: precomp should not have failed\n");
                fflush(stdout);
                flint_abort();
            }
        }
    }

cleanup:

    for (i = 0; i < r; i++)
    {
        fmpq_poly_clear(cQ + i);
        fmpq_poly_clear(bQ + i);
        fmpq_poly_clear(prod_bQ + i);
        fmpq_poly_clear(inv_prod_bQ + i);
    }
    flint_free(cQ);
    flint_free(bQ);
    flint_free(prod_bQ);
    flint_free(inv_prod_bQ);

    fmpz_poly_clear(a);
    fmpz_poly_clear(t);
    fmpz_poly_clear(t1);
    fmpq_poly_clear(aQ);
    fmpq_poly_clear(pQ);
    fmpq_poly_clear(G);
    fmpq_poly_clear(S);
}

TEST_FUNCTION_START(fmpz_poly_pfrac, state)
{
    slong i, j, k, tmul = 10;

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_poly_pfrac_t I;

        fmpz_poly_pfrac_init(I);

        for (j = 0; j < tmul; j++)
        {
            fmpz_poly_struct * b, * c;
            slong n = 2 + n_randint(state, 5);

            b = FLINT_ARRAY_ALLOC(n, fmpz_poly_struct);
            c = FLINT_ARRAY_ALLOC(n, fmpz_poly_struct);

            for (k = 0; k < n; k++)
            {
                fmpz_poly_init(c + k);
                fmpz_poly_init(b + k);
                do {
                    fmpz_poly_randtest(b + k, state, 2 + n_randint(state, 5),
                                                    2 + n_randint(state, 350));
                } while (b[k].length < 2);
            }

            _test_pfrac(c, I, b, n, state);

            for (k = 0; k < n; k++)
            {
                fmpz_poly_clear(c + k);
                fmpz_poly_clear(b + k);
            }
            flint_free(c);
            flint_free(b);
        }

        fmpz_poly_pfrac_clear(I);
    }

    TEST_FUNCTION_END(state);
}
