/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "acb.h"

TEST_FUNCTION_START(acb_lindep, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong len, k, i, j, prec, num, rank;
        acb_ptr vec;
        fmpz_mat_t rel, M;
        fmpz * coeffs;

        /* k random (small) entries, and len - k integer combinations of
           them: at least len - k relations with small coefficients */
        len = 1 + n_randint(state, 6);
        k = 1 + n_randint(state, len);
        prec = 64 + n_randint(state, 200);

        vec = _acb_vec_init(len);
        coeffs = _fmpz_vec_init(len);
        fmpz_mat_init(rel, len, len);
        fmpz_mat_init(M, len - k, len);

        for (i = 0; i < k; i++)
        {
            /* (generic entries: acb_randtest gives small dyadic numbers
               such as 1/2 and 1/16, which have further relations) */
            arb_urandom(acb_realref(vec + i), state, prec);
            arb_mul_2exp_si(acb_realref(vec + i), acb_realref(vec + i), (slong) n_randint(state, 7) - 3);
            if (n_randint(state, 2))
                arb_zero(acb_imagref(vec + i));
            else
                arb_urandom(acb_imagref(vec + i), state, prec);
            mag_zero(arb_radref(acb_realref(vec + i)));
            mag_zero(arb_radref(acb_imagref(vec + i)));
        }

        for (i = k; i < len; i++)
        {
            acb_zero(vec + i);
            for (j = 0; j < k; j++)
            {
                fmpz_randtest(coeffs + j, state, 1 + n_randint(state, 8));
                acb_addmul_fmpz(vec + i, vec + j, coeffs + j, prec);
                fmpz_set(fmpz_mat_entry(M, i - k, j), coeffs + j);
            }
            fmpz_set_si(fmpz_mat_entry(M, i - k, i), -1);
            if (n_randint(state, 2))
            {
                arb_add_error_2exp_si(acb_realref(vec + i), -prec + 2);
                arb_add_error_2exp_si(acb_imagref(vec + i), -prec + 2);
            }
        }

        num = acb_lindep(rel, vec, len, prec);

        /* every relation returned must hold numerically */
        for (i = 0; i < num; i++)
        {
            acb_t s;
            acb_init(s);
            /* (evaluated exactly, as acb_lindep validates when the
               entries span a small range of exponents; at precision prec
               the rounding error of a combination with large coefficients
               could hide a nonzero residual or not) */
            acb_zero(s);
            for (j = 0; j < len; j++)
                acb_addmul_fmpz(s, vec + j, fmpz_mat_entry(rel, i, j), ARF_PREC_EXACT);
            if (!acb_contains_zero(s) || _fmpz_vec_is_zero(fmpz_mat_row(rel, i), len))
            {
                flint_printf("FAIL: relation does not hold\n");
                _acb_vec_printd(vec, len, 20); flint_printf("\n");
                fmpz_mat_print_pretty(rel); flint_printf("\n");
                flint_abort();
            }
            acb_clear(s);
        }

        /* the planted relations should be recovered (the random entries
           are generic: no further relations with small coefficients) */
        rank = fmpz_mat_rank(M);
        if (num < rank)
        {
            flint_printf("FAIL: too few relations (num %wd, rank %wd)\n", num, rank);
            _acb_vec_printd(vec, len, 20); flint_printf("\n");
            fmpz_mat_print_pretty(M); flint_printf("\n");
            fmpz_mat_print_pretty(rel); flint_printf("\n");
            flint_abort();
        }

        _acb_vec_clear(vec, len);
        _fmpz_vec_clear(coeffs, len);
        fmpz_mat_clear(rel);
        fmpz_mat_clear(M);
    }

    TEST_FUNCTION_END(state);
}
