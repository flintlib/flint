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
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"

int
nmod_mat_is_in_howell_form(const nmod_mat_t A)
{
    slong *pivots;
    slong i, j, r;
    int numberpivots = 0;
    int prevrowzero = 0;
    mp_ptr extra_row;
    mp_limb_t g;

    if (nmod_mat_is_zero(A))
        return 1;

    pivots = flint_malloc(A->r * sizeof(slong));

    if (!nmod_mat_is_zero_row(A, 0))
    {
        for (j = 0; j < A->c; j++)
        {
            if (nmod_mat_entry(A, 0, j))
            {
                if ((A->mod).n % nmod_mat_entry(A, 0, j))
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
        if (!nmod_mat_is_zero_row(A, i))
        {
            if (prevrowzero)
            {
                flint_free(pivots);
                return 0;
            }
            for (j = 0; j < A->c; j++)
            {
                if (nmod_mat_entry(A, i, j))
                {
                    if (j <= pivots[numberpivots - 1])
                    {
                        flint_free(pivots);
                        return 0;
                    }
                    if ((A->mod).n % nmod_mat_entry(A, i, j))
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
            if (nmod_mat_entry(A, j, pivots[i]) >= nmod_mat_entry(A, i, pivots[i]))
            {
                flint_free(pivots);
                return 0;
            }
        }
    }
    extra_row = _nmod_vec_init(A->c);

    for (i = 0; i < numberpivots; i++)
    {
        g = n_gcd(A->mod.n, nmod_mat_entry(A, i, pivots[i]));

        if (g == 1)
        {
            continue;
        }

        g = A->mod.n/g;

        _nmod_vec_scalar_mul_nmod(extra_row, A->rows[i], A->c, g, A->mod);

        for ( j = pivots[i] + 1; j < A->c; j++)
        {
            if (extra_row[j])
            {
                for ( r = i; r < numberpivots; r++)
                {
                    if (pivots[r] == j)
                    {
                        if(!(extra_row[j] % nmod_mat_entry(A, r, pivots[r])))
                        {
                            g = extra_row[j]/nmod_mat_entry(A, r, pivots[r]);
                            _nmod_vec_scalar_addmul_nmod(extra_row, A->rows[r],
                                A->c, nmod_neg(g, A->mod), A->mod);
                        }

                    }
                }
            }
        }
        if (!_nmod_vec_is_zero(extra_row, A->c))
        {
            _nmod_vec_clear(extra_row);
            flint_free(pivots);
            return 0;
        }
    }
    _nmod_vec_clear(extra_row);
    flint_free(pivots);
    return 1;
}

TEST_FUNCTION_START(nmod_mat_howell_form, state)
{
    slong i;

    for (i = 0; i < 10000*flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, D;
        mp_limb_t mod;
        slong j, k, m, n, r1, r2;
        slong *perm;
        int equal;
        mp_limb_t c;

        mod = n_randtest_not_zero(state);

        m = n_randint(state, 20);
        do { n = n_randint(state, 20); } while (n > m);

        perm = _perm_init(2*m);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(D, 2*m, n, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_init_set(B, A);

        r1 = nmod_mat_howell_form(B);

        if (!nmod_mat_is_in_howell_form(B))
        {
            flint_printf("FAIL (malformed Howell form)\n");
            nmod_mat_print_pretty(A); flint_printf("\n\n");
            nmod_mat_print_pretty(B); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        _perm_randtest(perm, 2 * m, state);

        /* Concatenate the original matrix with the Howell form, scramble the rows,
           and check that the Howell form is the same */

        for (j = 0; j < m; j++)
        {
            do { c = n_randint(state, mod); } while ( n_gcd(c, mod) != 1);
            for (k = 0; k < n; k++)
                nmod_mat_entry(D, perm[j], k) =
                    nmod_mul(nmod_mat_entry(A, j, k), c, A->mod);
        }

        for (j = 0; j < m; j++)
        {
            do { c = n_randint(state, mod); } while ( n_gcd(c, mod) != 1);
            for (k = 0; k < n; k++)
                nmod_mat_entry(D, perm[m + j], k) =
                    nmod_mul(nmod_mat_entry(B, j, k), c, A->mod);
        }

        r2 = nmod_mat_howell_form(D);

        equal = (r1 == r2);

        if (equal)
        {
            for (j = 0; j < r1; j++)
                for (k = 0; k < n; k++)
                    equal = equal && (nmod_mat_entry(B, j, k) ==
                                        nmod_mat_entry(D, j, k));
            for (j = r1; j < 2*m; j++)
                for (k = 0; k < n; k++)
                    equal = equal && (nmod_mat_entry(D, j, k) == 0);
        }

        if (!equal)
        {
            flint_printf("FAIL (r1 = %wd, r2 = %wd)!\n", r1, r2);
            nmod_mat_print_pretty(A); flint_printf("\n\n");
            nmod_mat_print_pretty(B); flint_printf("\n\n");
            nmod_mat_print_pretty(D); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        _perm_clear(perm);

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(D);
    }

    TEST_FUNCTION_END(state);
}
