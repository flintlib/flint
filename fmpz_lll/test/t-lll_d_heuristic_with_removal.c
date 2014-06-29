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

    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

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
    fmpz_mat_t mat;
    fmpz_t bound;
    fmpz_lll_t fl;
    mp_bitcnt_t bits;

    FLINT_TEST_INIT(state);

    flint_printf("lll_d_heuristic_with_removal....");
    fflush(stdout);

    /* test using NTRU like matrices */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        ulong q;
        slong r, c;

        r = 2 * (n_randint(state, 50) + 1);
        c = r;

        fmpz_mat_init(mat, r, c);
        fmpz_init(bound);
        fmpz_lll_randtest(fl, state);

        bits = n_randint(state, 20) + 1;
        q = n_randint(state, 200) + 1;

        if (n_randint(state, 1))
            fmpz_mat_randntrulike(mat, state, bits, q);
        else
            fmpz_mat_randntrulike2(mat, state, bits, q);

        fmpz_randtest(bound, state, bits);
        if (fl->rt == GRAM)
        {
            fmpz_mat_gram(mat, mat);
            newd = fmpz_lll_d_heuristic_with_removal(mat, bound, fl);
            if (newd == -1)
            {
                fmpz_mat_clear(mat);
                fmpz_clear(bound);
                continue;
            }
            result =
                fmpz_mat_is_reduced_gram_with_removal(mat, fl->delta, fl->eta,
                                                      bound, newd);
        }
        else
        {
            newd = fmpz_lll_d_heuristic_with_removal(mat, bound, fl);
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
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong r, c;
        fmpz_mat_t gmat;

        r = n_randint(state, 20) + 1;
        c = r + 1;

        fmpz_mat_init(mat, r, c);
        fmpz_init(bound);
        fmpz_lll_randtest(fl, state);

        bits = n_randint(state, 200) + 1;

        fmpz_mat_randintrel(mat, state, bits);

        fmpz_randtest(bound, state, bits);
        if (fl->rt == GRAM)
        {
            fmpz_mat_init(gmat, r, r);
            fmpz_mat_gram(gmat, mat);
            newd = fmpz_lll_d_heuristic_with_removal(gmat, bound, fl);
            if (newd == -1)
            {
                fmpz_mat_clear(mat);
                fmpz_mat_clear(gmat);
                fmpz_clear(bound);
                continue;
            }
            result =
                fmpz_mat_is_reduced_gram_with_removal(gmat, fl->delta, fl->eta,
                                                      bound, newd);
        }
        else
        {
            newd = fmpz_lll_d_heuristic_with_removal(mat, bound, fl);
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
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong r, c;

        r = n_randint(state, 10) + 1;
        c = r;

        fmpz_mat_init(mat, r, c);
        fmpz_init(bound);
        fmpz_lll_randtest(fl, state);

        fmpz_mat_randajtai(mat, state, 0.5);

        fmpz_randtest(bound, state, n_randint(state, 200));
        if (fl->rt == GRAM)
        {
            fmpz_mat_gram(mat, mat);
            newd = fmpz_lll_d_heuristic_with_removal(mat, bound, fl);
            if (newd == -1)
            {
                fmpz_mat_clear(mat);
                fmpz_clear(bound);
                continue;
            }
            result =
                fmpz_mat_is_reduced_gram_with_removal(mat, fl->delta, fl->eta,
                                                      bound, newd);
        }
        else
        {
            newd = fmpz_lll_d_heuristic_with_removal(mat, bound, fl);
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
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        mp_bitcnt_t bits2;
        slong r, c;

        r = n_randint(state, 200) + 1;
        c = r;

        fmpz_mat_init(mat, r, c);
        fmpz_init(bound);
        fmpz_lll_randtest(fl, state);

        bits = n_randint(state, 200) + 1;
        bits2 = n_randint(state, 5) + 1;

        fmpz_mat_randsimdioph(mat, state, bits, bits2);

        fmpz_randtest(bound, state, bits);
        if (fl->rt == GRAM)
        {
            fmpz_mat_gram(mat, mat);
            newd = fmpz_lll_d_heuristic_with_removal(mat, bound, fl);
            if (newd == -1)
            {
                fmpz_mat_clear(mat);
                fmpz_clear(bound);
                continue;
            }
            result =
                fmpz_mat_is_reduced_gram_with_removal(mat, fl->delta, fl->eta,
                                                      bound, newd);
        }
        else
        {
            newd = fmpz_lll_d_heuristic_with_removal(mat, bound, fl);
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
