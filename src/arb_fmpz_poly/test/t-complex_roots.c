/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly_factor.h"
#include "fmpq_poly.h"
#include "acb.h"
#include "arb_poly.h"
#include "arith.h"
#include "arb_fmpz_poly.h"

void
check_roots(const fmpz_poly_t poly, acb_srcptr roots, slong prec)
{
    arb_ptr real;
    acb_ptr upper;
    arb_poly_t rpoly;
    arb_t lead;

    slong i, j, num_real, num_upper, deg;

    deg = fmpz_poly_degree(poly);

    num_real = 0;
    for (i = 0; i < deg; i++)
        if (acb_is_real(roots + i))
            num_real++;

    num_upper = (deg - num_real) / 2;

    real = _arb_vec_init(num_real);
    upper = _acb_vec_init(num_upper);
    arb_poly_init(rpoly);
    arb_init(lead);

    for (i = 0; i < num_real; i++)
        arb_set(real + i, acb_realref(roots + i));

    for (i = 0; i < num_upper; i++)
        acb_set(upper + i, roots + num_real + 2 * i);

    arb_poly_product_roots_complex(rpoly, real, num_real, upper, num_upper, prec);
    arb_set_fmpz(lead, poly->coeffs + deg);
    arb_poly_scalar_mul(rpoly, rpoly, lead, prec);

    if (!arb_poly_contains_fmpz_poly(rpoly, poly))
    {
        flint_printf("FAIL!\n");
        flint_printf("deg = %wd, num_real = %wd, num_upper = %wd\n\n", deg, num_real, num_upper);
        for (i = 0; i < deg; i++)
        {
            acb_printn(roots + i, 30, 0);
            flint_printf("\n");
        }

        flint_printf("\npoly:\n");
        fmpz_poly_print(poly); flint_printf("\n\n");

        flint_printf("rpoly:\n");
        arb_poly_printd(rpoly, 30); flint_printf("\n\n");
        flint_abort();
    }

    for (i = 0; i < deg; i++)
    {
        for (j = i + 1; j < deg; j++)
        {
            if (acb_overlaps(roots + i, roots + j))
            {
                flint_printf("FAIL! (isolation)\n");
                flint_printf("deg = %wd, num_real = %wd, num_upper = %wd\n\n", deg, num_real, num_upper);
                for (i = 0; i < deg; i++)
                {
                    acb_printn(roots + i, 30, 0);
                    flint_printf("\n");
                }
            }
        }
    }

    _arb_vec_clear(real, num_real);
    _acb_vec_clear(upper, num_upper);
    arb_poly_clear(rpoly);
    arb_clear(lead);
}

TEST_FUNCTION_START(arb_fmpz_poly_complex_roots, state)
{
    slong iter;

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t f, g;
        fmpq_poly_t h;
        fmpz_poly_factor_t fac;
        fmpz_t t;
        acb_ptr roots;
        slong i, j, n, deg, prec, num_factors;
        int flags;

        prec = 20 + n_randint(state, 1000);
        flags = 0; /* ARB_FMPZ_POLY_ROOTS_VERBOSE; */

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpq_poly_init(h);
        fmpz_init(t);
        fmpz_poly_one(f);

        num_factors = 1 + n_randint(state, 3);

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
                    arith_chebyshev_t_polynomial(g, n);
                    break;
                case 2:
                    arith_chebyshev_u_polynomial(g, n);
                    break;
                case 3:
                    arith_legendre_polynomial(h, n);
                    fmpq_poly_get_numerator(g, h);
                    break;
                case 4:
                    arith_cyclotomic_polynomial(g, n);
                    break;
                case 5:
                    arith_swinnerton_dyer_polynomial(g, n % 4);
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

        if (!fmpz_poly_is_zero(f))
        {
            fmpz_poly_factor_init(fac);
            fmpz_poly_factor_squarefree(fac, f);

            for (i = 0; i < fac->num; i++)
            {
                deg = fmpz_poly_degree(fac->p + i);
                roots = _acb_vec_init(deg);
                arb_fmpz_poly_complex_roots(roots, fac->p + i, flags, prec);
                check_roots(fac->p + i, roots, prec);
                _acb_vec_clear(roots, deg);
            }

            fmpz_poly_factor_clear(fac);
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpq_poly_clear(h);
        fmpz_clear(t);
    }

    /* divide by zero in D-K iteration */
    {
        fmpz_poly_t f;
        acb_ptr roots;

        fmpz_poly_init(f);
        roots = _acb_vec_init(6);

        fmpz_poly_set_str(f, "7  25402042698578724632715842150384812072165708201873598515430741927709434502441548044497536658582884343075517644457858427539077202038196334089177500850739920720211950909770464886080670804176754154963931850468594331583498 223893180314223240984084491888407141086932132045832889810633599526603312517647486257597642989561630621890634231277838275174064496550528961859898776882676932807021286 795640482736635886180831008959046125790483017543686711742925168363012610075751154572719364427044897243091617811844423305126322668210906163297113908677853824166006037 4675174622252362090426964349867694514717613116237971666947773542341728151885745558219487022927593367256013888828 8306993067199362953716207209584242447525211917662146954370874791469552420765132896500330153494280277905784766392 24405899409125496050877753569127928027275807342290477559552 28910098348753799787840354537641329148612896619499101692036");

        arb_fmpz_poly_complex_roots(roots, f, 0, 128);
        check_roots(f, roots, 512);

        _acb_vec_clear(roots, 6);
        fmpz_poly_clear(f);
    }

    TEST_FUNCTION_END(state);
}
