/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "fmpq_poly.h"
#include "arb.h"
#include "arb_poly.h"
#include "arith.h"
#include "arb_fmpz_poly.h"

void
check_real_roots(const fmpz_poly_t poly, arb_srcptr roots, slong num_real, slong prec)
{
    arb_t t;
    slong i, deg, num_actual;

    deg = fmpz_poly_degree(poly);

    num_actual = fmpz_poly_num_real_roots(poly);

    if (num_real != num_actual)
    {
        flint_printf("FAIL!\n");
        flint_printf("deg = %wd, num_real = %wd, actual = %wd\n\n", deg, num_real, num_actual);
        flint_abort();
    }

    arb_init(t);

    for (i = 0; i < num_real; i++)
    {
        if (arb_rel_accuracy_bits(roots + i) < prec)
        {
            flint_printf("FAIL! accuracy %wd (expected %wd)\n", arb_rel_accuracy_bits(roots + i), prec);
            flint_printf("\npoly:\n");
            fmpz_poly_print(poly); flint_printf("\n\n");
            arb_printd(roots + i, 30); flint_printf("\n\n");
            flint_abort();
        }

        arb_fmpz_poly_evaluate_arb(t, poly, roots + i, prec);

        if (!arb_contains_zero(t))
        {
            flint_printf("FAIL! f(root) != 0\n");
            flint_printf("\npoly:\n");
            fmpz_poly_print(poly); flint_printf("\n\n");
            arb_printd(roots + i, 30); flint_printf("\n\n");
            arb_printd(t, 30); flint_printf("\n\n");
            flint_abort();
        }
    }

    for (i = 0; i < num_real - 1; i++)
    {
        if (!arb_lt(roots + i, roots + i + 1))
        {
            flint_printf("FAIL! roots not isolated or ordered\n");

            for (i = 0; i < deg; i++)
            {
                arb_printd(roots + i, 30); flint_printf("\n\n");
                flint_printf("\n");
                flint_abort();
            }
        }
    }

    arb_clear(t);
}

static void fmpz_poly_squarefree_part(fmpz_poly_t res, fmpz_poly_t poly)
{

    if (poly->length == 0)
        fmpz_poly_zero(res);
    else if (poly->length == 1)
    {
        fmpz_poly_one(res);
    }
    else
    {
        fmpz_poly_t der, gcd;

        fmpz_poly_init(der);
        fmpz_poly_init(gcd);

        fmpz_poly_derivative(der, poly);
        fmpz_poly_gcd(gcd, poly, der);
        fmpz_poly_divexact(res, poly, gcd);

        fmpz_poly_clear(der);
        fmpz_poly_clear(gcd);

        if (res->length && fmpz_cmp_ui(res->coeffs + res->length - 1, 0) < 0)
        {
            fmpz_poly_neg(res, res);
        }
    }
}

TEST_FUNCTION_START(arb_fmpz_poly_real_roots, state)
{
    slong iter;

    static const char * regressions[3] = {
        "73  -1 -1 0 -1 -1 1 -1 0 0 0 0 0 1 0 0 0 1 -1 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 1 -1 -1 0 0 0 0 1 0 1 0 0 0 1 0 1 0 1 0 -1 -1 0 -1 0 0 -1 -1 0 1 0 0 0 0 0 0 -1 1",
        "5  3 0 -1 0 -2",
        "3  0 -5 -1",
    };

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t f, g;
        fmpq_poly_t h;
        fmpz_t t;
        arb_ptr roots;
        slong i, j, n, deg, prec, num_factors, num;
        int flags;

        prec = 20 + n_randint(state, 1000);
        flags = 0; /* ARB_FMPZ_POLY_ROOTS_VERBOSE; */

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpq_poly_init(h);
        fmpz_init(t);
        fmpz_poly_one(f);

        num_factors = 1 + n_randint(state, 3);

        if (iter < 3)
        {
            fmpz_poly_set_str(f, regressions[iter]);
        }
        else
        {
            for (i = 0; i < num_factors; i++)
            {
                n = n_randint(state, 18);

                switch (n_randint(state, 12))
                {
                    case 0:
                        fmpz_poly_zero(g);
                        for (j = 0; j <= n; j++)
                            fmpz_poly_set_coeff_ui(g, j, j+1);
                        break;
                    case 1:
                        fmpz_poly_chebyshev_t(g, n);
                        break;
                    case 2:
                        fmpz_poly_chebyshev_u(g, n);
                        break;
                    case 3:
                        fmpq_poly_legendre_p(h, n);
                        fmpq_poly_get_numerator(g, h);
                        break;
                    case 4:
                        fmpz_poly_cyclotomic(g, n);
                        break;
                    case 5:
                        fmpz_poly_swinnerton_dyer(g, n % 4);
                        break;
                    case 6:
                        arith_bernoulli_polynomial(h, n);
                        fmpq_poly_get_numerator(g, h);
                        break;
                    case 7:
                        fmpz_poly_zero(g);
                        fmpz_poly_fit_length(g, n+2);
                        arith_stirling_number_1_vec(g->coeffs, n+1, n+2);
                        _fmpz_poly_set_length(g, n+2);
                        fmpz_poly_shift_right(g, g, 1);
                        break;
                    case 8:
                        fmpq_poly_zero(h);
                        fmpq_poly_set_coeff_si(h, 0, 0);
                        fmpq_poly_set_coeff_si(h, 1, 1);
                        fmpq_poly_exp_series(h, h, n + 1);
                        fmpq_poly_get_numerator(g, h);
                        break;
                    case 9:
                        fmpz_poly_zero(g);
                        fmpz_poly_set_coeff_ui(g, 0, 1);
                        fmpz_poly_set_coeff_ui(g, 1, 100);
                        fmpz_poly_pow(g, g, n_randint(state, 5));
                        fmpz_poly_set_coeff_ui(g, n, 1);
                        break;
                    default:
                        fmpz_poly_randtest(g, state, 1 + n, 1 + n_randint(state, 300));
                        break;
                }

                fmpz_poly_mul(f, f, g);
            }
        }

        if (!fmpz_poly_is_zero(f))
        {
            fmpz_poly_squarefree_part(g, f);

            deg = fmpz_poly_degree(g);
            roots = _arb_vec_init(deg);

            num = arb_fmpz_poly_real_roots(roots, g, flags, prec);
            check_real_roots(g, roots, num, prec);

            _arb_vec_clear(roots, deg);
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpq_poly_clear(h);
        fmpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
