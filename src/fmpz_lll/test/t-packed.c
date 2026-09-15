/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "long_extras.h"
#include "fmpz_lll.h"

TEST_FUNCTION_START(fmpz_lll_packed, state)
{
    slong iter;

    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t B, C;
        fmpz_lll_packed_t P;
        slong d, n, bits, m, k, i, j, len;
        fmpz_t x, dot1, dot2;
        double * v1, * v2;
        slong e1, e2;

        d = 1 + n_randint(state, 8);
        n = 1 + n_randint(state, 12);
        bits = 1 + n_randint(state, 300);

        fmpz_mat_init(B, d, n);
        fmpz_mat_init(C, d, n);
        fmpz_init(x);
        fmpz_init(dot1);
        fmpz_init(dot2);
        fmpz_mat_randtest(B, state, bits);

        m = fmpz_lll_packed_limbs(B);
        if (n_randint(state, 2))
            m += n_randint(state, 3);

        fmpz_lll_packed_init(P, d, n, m);
        fmpz_lll_packed_set_fmpz_mat(P, B);

        /* random sequence of row operations, mirrored on B */
        for (k = 0; k < 30; k++)
        {
            int op = n_randint(state, 7);
            i = n_randint(state, d);
            j = n_randint(state, d);

            if (op == 0 && i != j)
            {
                fmpz_lll_packed_row_sub(P, i, j);
                _fmpz_vec_sub(fmpz_mat_row(B, i), fmpz_mat_row(B, i), fmpz_mat_row(B, j), n);
            }
            else if (op == 1 && i != j)
            {
                fmpz_lll_packed_row_add(P, i, j);
                _fmpz_vec_add(fmpz_mat_row(B, i), fmpz_mat_row(B, i), fmpz_mat_row(B, j), n);
            }
            else if (op == 2 && i != j)
            {
                slong xs = z_randtest(state);
                fmpz_lll_packed_row_submul_si(P, i, j, xs);
                _fmpz_vec_scalar_submul_si(fmpz_mat_row(B, i), fmpz_mat_row(B, j), n, xs);
            }
            else if (op == 3 && i != j)
            {
                ulong e = n_randint(state, 200);
                slong xs = z_randtest(state);
                fmpz_lll_packed_row_submul_si_2exp(P, i, j, xs, e);
                _fmpz_vec_scalar_submul_si_2exp(fmpz_mat_row(B, i), fmpz_mat_row(B, j), n, xs, e);
            }
            else if (op == 4 && i != j)
            {
                fmpz_randtest(x, state, 1 + n_randint(state, 200));
                fmpz_lll_packed_row_submul_fmpz(P, i, j, x);
                _fmpz_vec_scalar_submul_fmpz(fmpz_mat_row(B, i), fmpz_mat_row(B, j), n, x);
            }
            else if (op == 5)
            {
                fmpz_lll_packed_move_row(P, i, j);
                if (i > j)
                {
                    for (k = i; k > j; k--) fmpz_mat_swap_rows(B, NULL, k, k - 1);
                }
                else
                {
                    for (k = i; k < j; k++) fmpz_mat_swap_rows(B, NULL, k, k + 1);
                }
                k = 0;
            }
            else
            {
                fmpz_lll_packed_tighten(P, i);
                fmpz_lll_packed_maybe_shrink(P);
            }

            /* row bounds must be valid upper bounds */
            for (i = 0; i < d; i++)
            {
                slong b = 0;
                for (j = 0; j < n; j++)
                    b = FLINT_MAX(b, (slong) fmpz_bits(fmpz_mat_entry(B, i, j)));
                if (b > P->bits[i])
                    TEST_FUNCTION_FAIL("row bound not an upper bound\n");
            }
        }

        fmpz_lll_packed_get_fmpz_mat(C, P);

        if (!fmpz_mat_equal(B, C))
        {
            TEST_FUNCTION_FAIL("row operations differ\nB = %{fmpz_mat}\nC = %{fmpz_mat}\n", B, C);
        }

        /* dot product and double conversion */
        i = n_randint(state, d);
        j = n_randint(state, d);
        len = n_randint(state, n + 1);
        fmpz_lll_packed_dot(dot1, P, i, j, len);
        _fmpz_vec_dot(dot2, fmpz_mat_row(B, i), fmpz_mat_row(B, j), len);

        if (!fmpz_equal(dot1, dot2))
        {
            TEST_FUNCTION_FAIL("dot products differ: %{fmpz} vs %{fmpz}\n", dot1, dot2);
        }

        v1 = flint_malloc(sizeof(double) * n);
        v2 = flint_malloc(sizeof(double) * n);
        e1 = fmpz_lll_packed_get_d_vec_2exp(v1, P, i);
        e2 = _fmpz_vec_get_d_vec_2exp(v2, fmpz_mat_row(B, i), n);

        if (e1 != e2 || memcmp(v1, v2, sizeof(double) * n) != 0)
        {
            TEST_FUNCTION_FAIL("double conversion differs\n");
        }

        flint_free(v1);
        flint_free(v2);
        fmpz_lll_packed_clear(P);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_clear(x);
        fmpz_clear(dot1);
        fmpz_clear(dot2);
    }

    TEST_FUNCTION_END(state);
}
