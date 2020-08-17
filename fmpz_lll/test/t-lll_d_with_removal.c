/*
    Copyright (C) 2014 Abhinav Baid

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
#include "fmpz.h"
#include "fmpz_lll.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result = 1, newd;
    fmpz_mat_t mat, mat2, U;
    fmpz_t bound;
    fmpz_lll_t fl;
    flint_bitcnt_t bits;

    FLINT_TEST_INIT(state);

    flint_printf("lll_d_with_removal....");
    fflush(stdout);

    /* test using NTRU like matrices */
    for (i = 0; i < 1 * flint_test_multiplier(); i++)
    {
        ulong q;
        slong r, c;

        r = 2 * (n_randint(state, 25) + 1);
        c = r;

        fmpz_mat_init(mat, r, c);
        fmpz_init_set_ui(bound, 2);
        fmpz_lll_randtest(fl, state);

        bits = n_randint(state, 20) + 1;
        q = n_randint(state, 200) + 1;

        if (n_randint(state, 2))
            fmpz_mat_randntrulike(mat, state, bits, q);
        else
            fmpz_mat_randntrulike2(mat, state, bits, q);
        fmpz_mul_2exp(bound, bound, bits);

        if (fl->rt == GRAM)
        {
            if (n_randint(state, 2))
            {
                fmpz_mat_init(U, r, r);
                fmpz_mat_one(U);
                fmpz_mat_init(mat2, r, c);
                fmpz_mat_set(mat2, mat);
                fmpz_mat_gram(mat, mat);
                newd = fmpz_lll_d_with_removal(mat, U, bound, fl);
                fmpz_mat_mul(mat2, U, mat2);
                fmpz_mat_gram(mat2, mat2);
                result = fmpz_mat_equal(mat, mat2);
                if (!result)
                {
                    flint_printf
                        ("FAIL (randntrulike): gram matrices not equal!\n");
                    fmpz_mat_print_pretty(mat);
                    fmpz_mat_print_pretty(mat2);
                    abort();
                }
                fmpz_mat_clear(U);
                fmpz_mat_clear(mat2);
            }
            else
            {
                fmpz_mat_gram(mat, mat);
                newd = fmpz_lll_d_with_removal(mat, NULL, bound, fl);
            }
            result =
                fmpz_mat_is_reduced_gram_with_removal(mat, fl->delta, fl->eta,
                                                      bound, newd);
        }
        else
        {
            if (n_randint(state, 2))
            {
                fmpz_mat_init(U, r, r);
                fmpz_mat_one(U);
                fmpz_mat_init(mat2, r, c);
                fmpz_mat_set(mat2, mat);
                newd = fmpz_lll_d_with_removal(mat, U, bound, fl);
                fmpz_mat_mul(mat2, U, mat2);
                result = fmpz_mat_equal(mat, mat2);
                if (!result)
                {
                    flint_printf
                        ("FAIL (randntrulike): basis matrices not equal!\n");
                    fmpz_mat_print_pretty(mat);
                    fmpz_mat_print_pretty(mat2);
                    abort();
                }
                fmpz_mat_clear(U);
                fmpz_mat_clear(mat2);
            }
            else
            {
                newd = fmpz_lll_d_with_removal(mat, NULL, bound, fl);
            }
            result =
                fmpz_mat_is_reduced_with_removal(mat, fl->delta, fl->eta,
                                                 bound, newd);
        }

        if (!result)
        {
            flint_printf("FAIL (randntrulike):\n");
            fmpz_mat_print_pretty(mat);
            flint_printf("bits = %ld, i = %ld\n", bits, i);
            flint_printf("delta = %g, eta = %g\n", fl->delta, fl->eta);
            flint_printf("rep_type = %d\n", fl->rt);
            flint_printf("gram_type = %d\n", fl->gt);
            abort();
        }

        fmpz_mat_clear(mat);
        fmpz_clear(bound);
    }

    /* test using integer relations matrices */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        slong r, c;
        fmpz_mat_t gmat, gmat2;

        r = n_randint(state, 20) + 1;
        c = r + 1;

        fmpz_mat_init(mat, r, c);
        fmpz_init_set_ui(bound, 2);
        fmpz_lll_randtest(fl, state);

        bits = n_randint(state, 200) + 1;

        fmpz_mat_randintrel(mat, state, bits);
        fmpz_mul_2exp(bound, bound, bits);

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
                newd = fmpz_lll_d_with_removal(gmat, U, bound, fl);
                fmpz_mat_mul(mat2, U, mat);
                fmpz_mat_gram(gmat2, mat2);
                result = fmpz_mat_equal(gmat, gmat2);
                if (!result)
                {
                    flint_printf
                        ("FAIL (randintrel): gram matrices not equal!\n");
                    fmpz_mat_print_pretty(gmat);
                    fmpz_mat_print_pretty(gmat2);
                    abort();
                }
                fmpz_mat_clear(U);
                fmpz_mat_clear(mat2);
                fmpz_mat_clear(gmat2);
            }
            else
            {
                newd = fmpz_lll_d_with_removal(gmat, NULL, bound, fl);
            }
            result =
                fmpz_mat_is_reduced_gram_with_removal(gmat, fl->delta, fl->eta,
                                                      bound, newd);
        }
        else
        {
            if (n_randint(state, 2))
            {
                fmpz_mat_init(U, r, r);
                fmpz_mat_one(U);
                fmpz_mat_init(mat2, r, c);
                fmpz_mat_set(mat2, mat);
                newd = fmpz_lll_d_with_removal(mat, U, bound, fl);
                fmpz_mat_mul(mat2, U, mat2);
                result = fmpz_mat_equal(mat, mat2);
                if (!result)
                {
                    flint_printf
                        ("FAIL (randintrel): basis matrices not equal!\n");
                    fmpz_mat_print_pretty(mat);
                    fmpz_mat_print_pretty(mat2);
                    abort();
                }
                fmpz_mat_clear(U);
                fmpz_mat_clear(mat2);
            }
            else
            {
                newd = fmpz_lll_d_with_removal(mat, NULL, bound, fl);
            }
            result =
                fmpz_mat_is_reduced_with_removal(mat, fl->delta, fl->eta,
                                                 bound, newd);
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
            abort();
        }

        fmpz_mat_clear(mat);
        if (fl->rt == GRAM)
        {
            fmpz_mat_clear(gmat);
        }
        fmpz_clear(bound);
    }

    /* test using ajtai matrices */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        slong r, c;

        r = n_randint(state, 10) + 1;
        c = r;

        fmpz_mat_init(mat, r, c);
        fmpz_init_set_ui(bound, 4 * r);
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
                newd = fmpz_lll_d_with_removal(mat, U, bound, fl);
                fmpz_mat_mul(mat2, U, mat2);
                fmpz_mat_gram(mat2, mat2);
                result = fmpz_mat_equal(mat, mat2);
                if (!result)
                {
                    flint_printf
                        ("FAIL (randajtai): gram matrices not equal!\n");
                    fmpz_mat_print_pretty(mat);
                    fmpz_mat_print_pretty(mat2);
                    abort();
                }
                fmpz_mat_clear(U);
                fmpz_mat_clear(mat2);
            }
            else
            {
                fmpz_mat_gram(mat, mat);
                newd = fmpz_lll_d_with_removal(mat, NULL, bound, fl);
            }
            result =
                fmpz_mat_is_reduced_gram_with_removal(mat, fl->delta, fl->eta,
                                                      bound, newd);
        }
        else
        {
            if (n_randint(state, 2))
            {
                fmpz_mat_init(U, r, r);
                fmpz_mat_one(U);
                fmpz_mat_init(mat2, r, c);
                fmpz_mat_set(mat2, mat);
                newd = fmpz_lll_d_with_removal(mat, U, bound, fl);
                fmpz_mat_mul(mat2, U, mat2);
                result = fmpz_mat_equal(mat, mat2);
                if (!result)
                {
                    flint_printf
                        ("FAIL (randajtai): basis matrices not equal!\n");
                    fmpz_mat_print_pretty(mat);
                    fmpz_mat_print_pretty(mat2);
                    abort();
                }
                fmpz_mat_clear(U);
                fmpz_mat_clear(mat2);
            }
            else
            {
                newd = fmpz_lll_d_with_removal(mat, NULL, bound, fl);
            }
            result =
                fmpz_mat_is_reduced_with_removal(mat, fl->delta, fl->eta,
                                                 bound, newd);
        }

        if (!result)
        {
            flint_printf("FAIL (randajtai):\n");
            fmpz_mat_print_pretty(mat);
            flint_printf("i = %ld\n", i);
            flint_printf("delta = %g, eta = %g\n", fl->delta, fl->eta);
            flint_printf("rep_type = %d\n", fl->rt);
            flint_printf("gram_type = %d\n", fl->gt);
            abort();
        }

        fmpz_mat_clear(mat);
        fmpz_clear(bound);
    }

    /* test using simultaneous diophantine matrices */
    for (i = 0; i < 1 * flint_test_multiplier(); i++)
    {
        flint_bitcnt_t bits2;
        slong r, c;

        r = n_randint(state, 50) + 1;
        c = r;

        fmpz_mat_init(mat, r, c);
        fmpz_init_set_ui(bound, 2);
        fmpz_lll_randtest(fl, state);

        bits = n_randint(state, 200) + 1;
        bits2 = n_randint(state, 5) + 1;

        fmpz_mat_randsimdioph(mat, state, bits, bits2);
        fmpz_mul_2exp(bound, bound, bits2);

        if (fl->rt == GRAM)
        {
            if (n_randint(state, 2))
            {
                fmpz_mat_init(U, r, r);
                fmpz_mat_one(U);
                fmpz_mat_init(mat2, r, c);
                fmpz_mat_set(mat2, mat);
                fmpz_mat_gram(mat, mat);
                newd = fmpz_lll_d_with_removal(mat, U, bound, fl);
                fmpz_mat_mul(mat2, U, mat2);
                fmpz_mat_gram(mat2, mat2);
                result = fmpz_mat_equal(mat, mat2);
                if (!result)
                {
                    flint_printf
                        ("FAIL (randsimdioph): gram matrices not equal!\n");
                    fmpz_mat_print_pretty(mat);
                    fmpz_mat_print_pretty(mat2);
                    abort();
                }
                fmpz_mat_clear(U);
                fmpz_mat_clear(mat2);
            }
            else
            {
                fmpz_mat_gram(mat, mat);
                newd = fmpz_lll_d_with_removal(mat, NULL, bound, fl);
            }
            result =
                fmpz_mat_is_reduced_gram_with_removal(mat, fl->delta, fl->eta,
                                                      bound, newd);
        }
        else
        {
            if (n_randint(state, 2))
            {
                fmpz_mat_init(U, r, r);
                fmpz_mat_one(U);
                fmpz_mat_init(mat2, r, c);
                fmpz_mat_set(mat2, mat);
                newd = fmpz_lll_d_with_removal(mat, U, bound, fl);
                fmpz_mat_mul(mat2, U, mat2);
                result = fmpz_mat_equal(mat, mat2);
                if (!result)
                {
                    flint_printf
                        ("FAIL (randsimdioph): basis matrices not equal!\n");
                    fmpz_mat_print_pretty(mat);
                    fmpz_mat_print_pretty(mat2);
                    abort();
                }
                fmpz_mat_clear(U);
                fmpz_mat_clear(mat2);
            }
            else
            {
                newd = fmpz_lll_d_with_removal(mat, NULL, bound, fl);
            }
            result =
                fmpz_mat_is_reduced_with_removal(mat, fl->delta, fl->eta,
                                                 bound, newd);
        }

        if (!result)
        {
            flint_printf("FAIL (randsimdioph):\n");
            fmpz_mat_print_pretty(mat);
            flint_printf("bits = %ld, i = %ld\n", bits, i);
            flint_printf("delta = %g, eta = %g\n", fl->delta, fl->eta);
            flint_printf("rep_type = %d\n", fl->rt);
            flint_printf("gram_type = %d\n", fl->gt);
            abort();
        }

        fmpz_mat_clear(mat);
        fmpz_clear(bound);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
