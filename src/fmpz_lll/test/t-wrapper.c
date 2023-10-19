/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"
#include "fmpz_lll.h"

TEST_FUNCTION_START(fmpz_lll_wrapper, state)
{
    int i, result = 1;
    fmpz_mat_t mat, mat2, U;
    fmpz_lll_t fl;
    flint_bitcnt_t bits;

    /* test using NTRU like matrices */
    for (i = 0; i < 1 * flint_test_multiplier(); i++)
    {
        ulong q;
        slong r, c;

        r = 2 * (n_randint(state, 25) + 1);
        c = r;

        fmpz_mat_init(mat, r, c);
        fmpz_lll_randtest(fl, state);

        bits = n_randint(state, 20) + 1;
        q = n_randint(state, 200) + 1;

        if (n_randint(state, 2))
            fmpz_mat_randntrulike(mat, state, bits, q);
        else
            fmpz_mat_randntrulike2(mat, state, bits, q);

        if (fl->rt == GRAM)
        {
            if (n_randint(state, 2))
            {
                fmpz_mat_init(U, r, r);
                fmpz_mat_one(U);
                fmpz_mat_init(mat2, r, c);
                fmpz_mat_set(mat2, mat);
                fmpz_mat_gram(mat, mat);
                fmpz_lll_wrapper(mat, U, fl);
                fmpz_mat_mul(mat2, U, mat2);
                fmpz_mat_gram(mat2, mat2);
                result = fmpz_mat_equal(mat, mat2);
                if (!result)
                {
                    flint_printf
                        ("FAIL (randntrulike): gram matrices not equal!\n");
                    fmpz_mat_print_pretty(mat);
                    fmpz_mat_print_pretty(mat2);
                    fflush(stdout);
                    flint_abort();
                }
                fmpz_mat_clear(U);
                fmpz_mat_clear(mat2);
            }
            else
            {
                fmpz_mat_gram(mat, mat);
                fmpz_lll_wrapper(mat, NULL, fl);
            }
            result = fmpz_mat_is_reduced_gram(mat, fl->delta, fl->eta);
        }
        else
        {
            if (n_randint(state, 2))
            {
                fmpz_mat_init(U, r, r);
                fmpz_mat_one(U);
                fmpz_mat_init(mat2, r, c);
                fmpz_mat_set(mat2, mat);
                fmpz_lll_wrapper(mat, U, fl);
                fmpz_mat_mul(mat2, U, mat2);
                result = fmpz_mat_equal(mat, mat2);
                if (!result)
                {
                    flint_printf
                        ("FAIL (randntrulike): basis matrices not equal!\n");
                    fmpz_mat_print_pretty(mat);
                    fmpz_mat_print_pretty(mat2);
                    fflush(stdout);
                    flint_abort();
                }
                fmpz_mat_clear(U);
                fmpz_mat_clear(mat2);
            }
            else
            {
                fmpz_lll_wrapper(mat, NULL, fl);
            }
            result = fmpz_mat_is_reduced(mat, fl->delta, fl->eta);
        }

        if (!result)
        {
            flint_printf("FAIL (randntrulike):\n");
            fmpz_mat_print_pretty(mat);
            flint_printf("bits = %ld, i = %ld\n", bits, i);
            flint_printf("delta = %g, eta = %g\n", fl->delta, fl->eta);
            flint_printf("rep_type = %d\n", fl->rt);
            flint_printf("gram_type = %d\n", fl->gt);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(mat);
    }

    /* test using integer relations matrices */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        slong r, c;
        fmpz_mat_t gmat, gmat2;

        r = n_randint(state, 20) + 1;
        c = r + 1;

        fmpz_mat_init(mat, r, c);
        fmpz_lll_randtest(fl, state);

        bits = n_randint(state, 200) + 1;

        fmpz_mat_randintrel(mat, state, bits);

        if (fl->rt == GRAM)
        {
            fmpz_mat_init(gmat, r, r);
            fmpz_mat_gram(gmat, mat);
            if (n_randint(state, 2))
            {
                fmpz_mat_init(U, r, r);
                fmpz_mat_one(U);
                fmpz_mat_init(mat2, r, c);
                fmpz_mat_init(gmat2, r, r);
                fmpz_lll_wrapper(gmat, U, fl);
                fmpz_mat_mul(mat2, U, mat);
                fmpz_mat_gram(gmat2, mat2);
                result = fmpz_mat_equal(gmat, gmat2);
                if (!result)
                {
                    flint_printf
                        ("FAIL (randintrel): gram matrices not equal!\n");
                    fmpz_mat_print_pretty(gmat);
                    fmpz_mat_print_pretty(gmat2);
                    fflush(stdout);
                    flint_abort();
                }
                fmpz_mat_clear(U);
                fmpz_mat_clear(mat2);
                fmpz_mat_clear(gmat2);
            }
            else
            {
                fmpz_lll_wrapper(gmat, NULL, fl);
            }
            result = fmpz_mat_is_reduced_gram(gmat, fl->delta, fl->eta);
        }
        else
        {
            if (n_randint(state, 2))
            {
                fmpz_mat_init(U, r, r);
                fmpz_mat_one(U);
                fmpz_mat_init(mat2, r, c);
                fmpz_mat_set(mat2, mat);
                fmpz_lll_wrapper(mat, U, fl);
                fmpz_mat_mul(mat2, U, mat2);
                result = fmpz_mat_equal(mat, mat2);
                if (!result)
                {
                    flint_printf
                        ("FAIL (randintrel): basis matrices not equal!\n");
                    fmpz_mat_print_pretty(mat);
                    fmpz_mat_print_pretty(mat2);
                    fflush(stdout);
                    flint_abort();
                }
                fmpz_mat_clear(U);
                fmpz_mat_clear(mat2);
            }
            else
            {
                fmpz_lll_wrapper(mat, NULL, fl);
            }
            result = fmpz_mat_is_reduced(mat, fl->delta, fl->eta);
        }

        if (!result)
        {
            flint_printf("FAIL (randintrel):\n");
            fmpz_mat_print_pretty(mat);
            if (fl->rt == GRAM)
            {
                fmpz_mat_print_pretty(gmat);
            }
            flint_printf("bits = %ld, i = %ld\n", bits, i);
            flint_printf("delta = %g, eta = %g\n", fl->delta, fl->eta);
            flint_printf("rep_type = %d\n", fl->rt);
            flint_printf("gram_type = %d\n", fl->gt);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(mat);
        if (fl->rt == GRAM)
        {
            fmpz_mat_clear(gmat);
        }
    }

    /* test using ajtai matrices */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        slong r, c;

        r = n_randint(state, 10) + 1;
        c = r;

        fmpz_mat_init(mat, r, c);
        fmpz_lll_randtest(fl, state);

        fmpz_mat_randajtai(mat, state, 0.5);

        if (fl->rt == GRAM)
        {
            if (n_randint(state, 2))
            {
                fmpz_mat_init(U, r, r);
                fmpz_mat_one(U);
                fmpz_mat_init(mat2, r, c);
                fmpz_mat_set(mat2, mat);
                fmpz_mat_gram(mat, mat);
                fmpz_lll_wrapper(mat, U, fl);
                fmpz_mat_mul(mat2, U, mat2);
                fmpz_mat_gram(mat2, mat2);
                result = fmpz_mat_equal(mat, mat2);
                if (!result)
                {
                    flint_printf
                        ("FAIL (randajtai): gram matrices not equal!\n");
                    fmpz_mat_print_pretty(mat);
                    fmpz_mat_print_pretty(mat2);
                    fflush(stdout);
                    flint_abort();
                }
                fmpz_mat_clear(U);
                fmpz_mat_clear(mat2);
            }
            else
            {
                fmpz_mat_gram(mat, mat);
                fmpz_lll_wrapper(mat, NULL, fl);
            }
            result = fmpz_mat_is_reduced_gram(mat, fl->delta, fl->eta);
        }
        else
        {
            if (n_randint(state, 2))
            {
                fmpz_mat_init(U, r, r);
                fmpz_mat_one(U);
                fmpz_mat_init(mat2, r, c);
                fmpz_mat_set(mat2, mat);
                fmpz_lll_wrapper(mat, U, fl);
                fmpz_mat_mul(mat2, U, mat2);
                result = fmpz_mat_equal(mat, mat2);
                if (!result)
                {
                    flint_printf
                        ("FAIL (randajtai): basis matrices not equal!\n");
                    fmpz_mat_print_pretty(mat);
                    fmpz_mat_print_pretty(mat2);
                    fflush(stdout);
                    flint_abort();
                }
                fmpz_mat_clear(U);
                fmpz_mat_clear(mat2);
            }
            else
            {
                fmpz_lll_wrapper(mat, NULL, fl);
            }
            result = fmpz_mat_is_reduced(mat, fl->delta, fl->eta);
        }

        if (!result)
        {
            flint_printf("FAIL (randajtai):\n");
            fmpz_mat_print_pretty(mat);
            flint_printf("i = %ld\n", i);
            flint_printf("delta = %g, eta = %g\n", fl->delta, fl->eta);
            flint_printf("rep_type = %d\n", fl->rt);
            flint_printf("gram_type = %d\n", fl->gt);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(mat);
    }

    /* test using simultaneous diophantine matrices */
    for (i = 0; i < 1 * flint_test_multiplier(); i++)
    {
        flint_bitcnt_t bits2;
        slong r, c;

        r = n_randint(state, 50) + 1;
        c = r;

        fmpz_mat_init(mat, r, c);
        fmpz_lll_randtest(fl, state);

        bits = n_randint(state, 200) + 1;
        bits2 = n_randint(state, 5) + 1;

        fmpz_mat_randsimdioph(mat, state, bits, bits2);

        if (fl->rt == GRAM)
        {
            if (n_randint(state, 2))
            {
                fmpz_mat_init(U, r, r);
                fmpz_mat_one(U);
                fmpz_mat_init(mat2, r, c);
                fmpz_mat_set(mat2, mat);
                fmpz_mat_gram(mat, mat);
                fmpz_lll_wrapper(mat, U, fl);
                fmpz_mat_mul(mat2, U, mat2);
                fmpz_mat_gram(mat2, mat2);
                result = fmpz_mat_equal(mat, mat2);
                if (!result)
                {
                    flint_printf
                        ("FAIL (randsimdioph): gram matrices not equal!\n");
                    fmpz_mat_print_pretty(mat);
                    fmpz_mat_print_pretty(mat2);
                    fflush(stdout);
                    flint_abort();
                }
                fmpz_mat_clear(U);
                fmpz_mat_clear(mat2);
            }
            else
            {
                fmpz_mat_gram(mat, mat);
                fmpz_lll_wrapper(mat, NULL, fl);
            }
            result = fmpz_mat_is_reduced_gram(mat, fl->delta, fl->eta);
        }
        else
        {
            if (n_randint(state, 2))
            {
                fmpz_mat_init(U, r, r);
                fmpz_mat_one(U);
                fmpz_mat_init(mat2, r, c);
                fmpz_mat_set(mat2, mat);
                fmpz_lll_wrapper(mat, U, fl);
                fmpz_mat_mul(mat2, U, mat2);
                result = fmpz_mat_equal(mat, mat2);
                if (!result)
                {
                    flint_printf
                        ("FAIL (randsimdioph): basis matrices not equal!\n");
                    fmpz_mat_print_pretty(mat);
                    fmpz_mat_print_pretty(mat2);
                    fflush(stdout);
                    flint_abort();
                }
                fmpz_mat_clear(U);
                fmpz_mat_clear(mat2);
            }
            else
            {
                fmpz_lll_wrapper(mat, NULL, fl);
            }
            result = fmpz_mat_is_reduced(mat, fl->delta, fl->eta);
        }

        if (!result)
        {
            flint_printf("FAIL (randsimdioph):\n");
            fmpz_mat_print_pretty(mat);
            flint_printf("bits = %ld, i = %ld\n", bits, i);
            flint_printf("delta = %g, eta = %g\n", fl->delta, fl->eta);
            flint_printf("rep_type = %d\n", fl->rt);
            flint_printf("gram_type = %d\n", fl->gt);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(mat);
    }

    TEST_FUNCTION_END(state);
}
