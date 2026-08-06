/*
    Copyright (C) 2014 Abhinav Baid
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"
#include "fmpz_lll.h"

/* todo: test some examples where _d or _mpfr specifically fails */
TEST_FUNCTION_START(fmpz_lll_is_reduced, state)
{
    int i;
    fmpz_mat_t mat, gmat;
    fmpz_lll_t fl;
    flint_bitcnt_t bits;
    int reduced_d, reduced_mpfr, reduced;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_lll_randtest(fl, state);

        switch (n_randint(state, 4))
        {
            case 0:
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
                }
                break;
            case 1:
                {
                    slong r, c;

                    r = n_randint(state, 20) + 1;
                    c = r + 1;

                    fmpz_mat_init(mat, r, c);
                    bits = n_randint(state, 200) + 1;

                    fmpz_mat_randintrel(mat, state, bits);
                }
                break;
            case 2:
                {
                    slong r, c;

                    r = n_randint(state, 10) + 1;
                    c = r;

                    fmpz_mat_init(mat, r, c);
                    fmpz_lll_randtest(fl, state);

                    fmpz_mat_randajtai(mat, state, 0.5);
                }
                break;
            default:
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
                }
                break;
        }

        if (fl->rt == GRAM)
        {
            fmpz_mat_init(gmat, mat->r, mat->r);
            fmpz_mat_gram(gmat, mat);

            if (n_randint(state, 2))
                fmpz_lll_wrapper(gmat, NULL, fl);

            reduced = fmpz_mat_is_reduced_gram(gmat, fl->delta, fl->eta);
            reduced_d = fmpz_lll_is_reduced_d(gmat, fl);
            reduced_mpfr = fmpz_lll_is_reduced_mpfr(gmat, fl, 128);

            fmpz_mat_clear(gmat);
        }
        else
        {
            if (n_randint(state, 2))
                fmpz_lll_wrapper(mat, NULL, fl);

            reduced = fmpz_mat_is_reduced(mat, fl->delta, fl->eta);
            reduced_d = fmpz_lll_is_reduced_d(mat, fl);
            reduced_mpfr = fmpz_lll_is_reduced_mpfr(mat, fl, 128);
        }

        if (reduced_d && !reduced)
        {
            flint_printf("FAIL:\n");
            fmpz_mat_print_pretty(mat);
            flint_printf("delta = %g, eta = %g\n", fl->delta, fl->eta);
            flint_printf("rep_type = %d\n", fl->rt);
            flint_printf("gram_type = %d\n", fl->gt);
            flint_printf("reduced_d = %d, reduced = %d\n", reduced_d, reduced);
            fflush(stdout);
            flint_abort();
        }

        if (reduced_mpfr && !reduced)
        {
            flint_printf("FAIL:\n");
            fmpz_mat_print_pretty(mat);
            flint_printf("delta = %g, eta = %g\n", fl->delta, fl->eta);
            flint_printf("rep_type = %d\n", fl->rt);
            flint_printf("gram_type = %d\n", fl->gt);
            flint_printf("reduced_mpfr = %d, reduced = %d\n", reduced_mpfr, reduced);
            fflush(stdout);
            flint_abort();
        }


        fmpz_mat_clear(mat);
    }

    TEST_FUNCTION_END(state);
}
