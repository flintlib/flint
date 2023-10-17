/*
    Copyright (C) 2015 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "perm.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

int
fmpz_mat_mod_is_in_howell_form(const fmpz_mat_t A, const fmpz_t mod)
{
    slong *pivots;
    slong i, j, r;
    int numberpivots = 0;
    int prevrowzero = 0;
    fmpz * extra_row;
    fmpz_t g;

    if (fmpz_mat_is_empty(A))
        return 1;

    pivots = flint_malloc(A->r * sizeof(slong));

    if (!fmpz_mat_is_zero_row(A, 0))
    {
        for (j = 0; j < A->c; j++)
        {
            if (!fmpz_is_zero(fmpz_mat_entry(A, 0, j)))
            {
                if (!fmpz_divisible(mod, fmpz_mat_entry(A, 0, j)))
                {
                    flint_free(pivots);
                    return 0;
                }
                pivots[numberpivots] = j;
                numberpivots++;
                break;
            }
        }
    }
    else
    {
        prevrowzero = 1;
    }

    for (i = 1; i < A->r; i++)
    {
        if (!fmpz_mat_is_zero_row(A, i))
        {
            if (prevrowzero)
            {
                flint_free(pivots);
                return 0;
            }
            for (j = 0; j < A->c; j++)
            {
                if (!fmpz_is_zero(fmpz_mat_entry(A, i, j)))
                {
                    if (j <= pivots[numberpivots - 1])
                    {
                        flint_free(pivots);
                        return 0;
                    }
                    if (!fmpz_divisible(mod, fmpz_mat_entry(A, i, j)))
                    {
                        flint_free(pivots);
                        return 0;
                    }
                    pivots[numberpivots] = j;
                    numberpivots++;
                    j = A->c;
                }
            }
        }
        else
        {
            prevrowzero = 1;
        }
    }

    for (i = 1; i < numberpivots; i++)
    {
        for (j = 0; j < i; j++)
        {
            if (fmpz_cmp(fmpz_mat_entry(A, j, pivots[i]), fmpz_mat_entry(A, i, pivots[i])) >= 0)
            {
                flint_free(pivots);
                return 0;
            }
        }
    }
    extra_row = _fmpz_vec_init(A->c);

    fmpz_init(g);

    for (i = 0; i < numberpivots; i++)
    {
        fmpz_gcd(g, mod, fmpz_mat_entry(A, i, pivots[i]));
        fmpz_divexact(g, mod, g);
        _fmpz_vec_scalar_mul_fmpz(extra_row, A->rows[i], A->c, g);
        _fmpz_vec_scalar_mod_fmpz(extra_row, extra_row, A->c, mod);

        for ( j = pivots[i] + 1; j < A->c; j++)
        {
            if (!fmpz_is_zero(extra_row + j))
            {
                for ( r = i; r < numberpivots; r++)
                {
                    if (pivots[r] == j)
                    {
                        if (fmpz_divisible(extra_row + j, fmpz_mat_entry(A, r, pivots[r])))
                        {
                            fmpz_divexact(g, extra_row + j, fmpz_mat_entry(A, r, pivots[r]));
                            fmpz_neg(g, g);
                            _fmpz_vec_scalar_addmul_fmpz(extra_row, A->rows[r], A->c, g);
                        }
                    }
                }
            }
        }

        _fmpz_vec_scalar_mod_fmpz(extra_row, extra_row, A->c, mod);

        if (!_fmpz_vec_is_zero(extra_row, A->c))
        {
            _fmpz_vec_clear(extra_row, A->c);
            flint_free(pivots);
            fmpz_clear(g);
            return 0;
        }
    }
    _fmpz_vec_clear(extra_row, A->c);
    flint_free(pivots);
    fmpz_clear(g);
    return 1;
}

TEST_FUNCTION_START(fmpz_mat_howell_form_mod, state)
{
    slong i;

    for (i = 0; i < 10000*flint_test_multiplier(); i++)
    {
        fmpz_mat_t A, B, D;
        fmpz_t mod;
        fmpz_t t, c;
        slong j, k, m, n, r1, r2;
        slong *perm;
        int equal;

        fmpz_init(mod);
        fmpz_init(t);
        fmpz_init(c);

        do { fmpz_randtest_unsigned(mod, state, 10); } while (fmpz_is_zero(mod));

        m = n_randint(state, 20);
        do { n = n_randint(state, 20); } while (n > m);

        perm = _perm_init(2*m);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(D, 2*m, n);

        fmpz_mat_randtest(A, state, 10);
        fmpz_mat_init_set(B, A);

        r1 = fmpz_mat_howell_form_mod(B, mod);

        if (!fmpz_mat_mod_is_in_howell_form(B, mod))
        {
            flint_printf("FAIL (malformed Howell form)\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(B); flint_printf("\n\n");
            flint_printf("Modulus: ");
            fmpz_print(mod);
            flint_printf("\n\n");

            fflush(stdout);
            flint_abort();
        }

        _perm_randtest(perm, 2 * m, state);

        /* Concatenate the original matrix with the Howell form, scramble the rows,
           and check that the Howell form is the same */

        for (j = 0; j < m; j++)
        {
            while (1)
            {
                fmpz_randtest_mod(c, state, mod);
                fmpz_gcd(t, c, mod);
                if (fmpz_is_one(t)) { break; }
            }

            for (k = 0; k < n; k++)
                fmpz_mul(fmpz_mat_entry(D, perm[j], k), fmpz_mat_entry(A, j, k), c);
        }

        for (j = 0; j < m; j++)
        {
            while (1)
            {
                fmpz_randtest_mod(c, state, mod);
                fmpz_gcd(t, c, mod);
                if (fmpz_is_one(t)) { break; }
            }

            for (k = 0; k < n; k++)
                fmpz_mul(fmpz_mat_entry(D, perm[m + j], k), fmpz_mat_entry(B, j, k), c);
        }

        r2 = fmpz_mat_howell_form_mod(D, mod);

        equal = (r1 == r2);

        if (equal)
        {
            for (j = 0; j < r1; j++)
                for (k = 0; k < n; k++)
                    equal = equal && fmpz_equal(fmpz_mat_entry(B, j, k),fmpz_mat_entry(D, j, k));
            for (j = r1; j < 2*m; j++)
                for (k = 0; k < n; k++)
                    equal = equal && fmpz_is_zero(fmpz_mat_entry(D, j, k));
        }

        if (!equal)
        {
            flint_printf("FAIL (r1 = %wd, r2 = %wd)!\n", r1, r2);
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(B); flint_printf("\n\n");
            fmpz_mat_print_pretty(D); flint_printf("\n\n");
            flint_printf("Modulus: ");
            fmpz_print(mod);

            fflush(stdout);
            flint_abort();
        }

        _perm_clear(perm);

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(D);

        fmpz_clear(mod);
        fmpz_clear(t);
        fmpz_clear(c);
    }

    TEST_FUNCTION_END(state);
}
