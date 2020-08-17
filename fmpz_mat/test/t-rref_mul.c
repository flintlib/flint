/*
    Copyright (C) 2010-2012 Fredrik Johansson
    Copyright (C) 2014 Alex J. Best

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
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "perm.h"
#include "ulong_extras.h"

int
main(void)
{
    slong iter;
    FLINT_TEST_INIT(state);

    flint_printf("rref_mul....");
    fflush(stdout);

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, R, B, R2, R3;
        fmpz_t den, c, den2;
        slong j, k, m, n, b, d, r, rank1, rank2;
        slong * perm;
        int equal;

        m = 1 + n_randint(state, 10);
        n = 1 + n_randint(state, 10);
        r = n_randint(state, FLINT_MIN(m, n) + 1);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(R, m, n);
        fmpz_mat_init(R2, m, n);
        fmpz_mat_init(R3, 2 * m, n);
        fmpz_mat_init(B, 2 * m, n);

        fmpz_init(c);
        fmpz_init(den);
        fmpz_init(den2);

        perm = _perm_init(2 * m);

        /* sparse */
        b = 1 + n_randint(state, 10) * n_randint(state, 10);
        d = n_randint(state, 2*m*n + 1);
        fmpz_mat_randrank(A, state, r, b);

        /* dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, d);

        rank1 = fmpz_mat_rref_mul(R, den, A);

        if (r != rank1)
        {
            flint_printf("FAIL wrong rank! (r = %wd, rank1 = %wd)!\n", r, rank1);
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(R); flint_printf("\n\n");
            abort();
        }

        if (!fmpz_mat_is_in_rref_with_rank(R, den, rank1))
        {
            flint_printf("FAIL matrix not in rref!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(R); flint_printf("\n\n");
            abort();
        }

        /* check rref is the same as when computed classically */

        rank2 = fmpz_mat_rref_fflu(R2, den2, A);
        equal = (rank1 == rank2);

        if (equal)
        {
            fmpz_mat_scalar_mul_fmpz(R, R, den2);
            fmpz_mat_scalar_mul_fmpz(R2, R2, den);

            for (j = 0; j < m; j++)
                for (k = 0; k < n; k++)
                    equal = equal && fmpz_equal(fmpz_mat_entry(R, j, k),
                            fmpz_mat_entry(R2, j, k));

            if (!fmpz_is_zero(den2)) /* TODO should den2 be able to be? */
                fmpz_mat_scalar_divexact_fmpz(R, R, den2);
        }

        if (!equal)
        {
            flint_printf("FAIL different to classical (rank1 = %wd, rank2 = %wd)!\n", rank1, rank2);
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(R); flint_printf("\n\n");
            fmpz_mat_print_pretty(R2); flint_printf("\n\n");
            abort();
        }

        /* Concatenate the original matrix with the rref, scramble the rows,
            and check that the rref is the same */
        _perm_randtest(perm, 2 * m, state);

        for (j = 0; j < m; j++)
        {
            fmpz_randtest_not_zero(c, state, 5);
            for (k = 0; k < n; k++)
                fmpz_mul(fmpz_mat_entry(B, perm[j], k), fmpz_mat_entry(A, j, k), c);
        }

        for (j = 0; j < m; j++)
        {
            fmpz_randtest_not_zero(c, state, 5);
            for (k = 0; k < n; k++)
                fmpz_mul(fmpz_mat_entry(B, perm[m + j], k), fmpz_mat_entry(R, j, k), c);
        }

        rank2 = fmpz_mat_rref_mul(R3, den2, B);
        equal = (rank1 == rank2);

        if (equal)
        {
            fmpz_mat_scalar_mul_fmpz(R, R, den2);
            fmpz_mat_scalar_mul_fmpz(R3, R3, den);

            for (j = 0; j < rank2; j++)
                for (k = 0; k < n; k++)
                    equal = equal &&
                        fmpz_equal(fmpz_mat_entry(R, j, k), fmpz_mat_entry(R3, j, k));
            for (j = rank2; j < 2 * rank2; j++)
                for (k = 0; k < n; k++)
                    equal = equal && fmpz_is_zero(fmpz_mat_entry(R3, j, k));
        }

        if (!equal)
        {
            flint_printf("FAIL (rank1 = %wd, rank2 = %wd)!\n", rank1, rank2);
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(R); flint_printf("\n\n");
            fmpz_mat_print_pretty(R3); flint_printf("\n\n");
            abort();
        }

        fmpz_clear(c);
        fmpz_clear(den);
        fmpz_clear(den2);

        _perm_clear(perm);

        fmpz_mat_clear(B);
        fmpz_mat_clear(R3);
        fmpz_mat_clear(R2);
        fmpz_mat_clear(R);
        fmpz_mat_clear(A);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
