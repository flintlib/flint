/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"

/* Sets p to a random squarefree polynomial with a known number of real
   roots (returned), of which *n01 lie in (0, 1). The polynomial is chosen
   to exercise the different algorithms used by fmpz_poly_num_real_roots
   and fmpz_poly_num_real_roots_0_1. */
static slong
_randtest_known_roots(fmpz_poly_t p, slong * n01, flint_rand_t state)
{
    fmpz_poly_t q, lin;
    fmpz_t a, b;
    slong m, j, nr, bits, len;

    fmpz_poly_init(q);
    fmpz_poly_init(lin);
    fmpz_init(a);
    fmpz_init(b);

    do
    {
        nr = 0;
        *n01 = 0;

        switch (n_randint(state, 5))
        {
            case 0: /* products of linear factors times a factor with no real roots */
                m = n_randint(state, 30);
                bits = 1 + n_randint(state, n_randint(state, 2) ? 8 : 100);
                len = n_randint(state, 2) ? 1 + n_randint(state, 5) : 1 + n_randint(state, 40);

                fmpz_poly_randtest_no_real_root(p, state, len, 1 + n_randint(state, 200));
                if (fmpz_poly_is_zero(p))
                    fmpz_poly_one(p);

                for (j = 0; j < m; j++)
                {
                    fmpz_randtest_not_zero(a, state, bits);
                    fmpz_abs(a, a);
                    if (n_randint(state, 4) == 0)
                    {
                        /* roots at, close to, or inside [0, 1] */
                        fmpz_randm(b, state, a);
                        if (n_randint(state, 8) == 0)
                            fmpz_zero(b);
                        else if (n_randint(state, 8) == 0)
                            fmpz_set(b, a);
                    }
                    else
                        fmpz_randtest(b, state, bits);

                    fmpz_poly_set_coeff_fmpz(lin, 1, a);
                    fmpz_neg(b, b);
                    fmpz_poly_set_coeff_fmpz(lin, 0, b);
                    fmpz_neg(b, b);
                    fmpz_poly_mul(p, p, lin);
                    nr++;
                    /* root b/a with a > 0: in (0, 1) iff 0 < b < a */
                    if (fmpz_sgn(b) > 0 && fmpz_cmp(b, a) < 0)
                        (*n01)++;
                }
                break;

            case 1: /* Chebyshev T_m: m real roots in (-1, 1), scaled by 2^k */
                m = 1 + n_randint(state, 120);
                fmpz_poly_chebyshev_t(p, m);
                j = n_randint(state, 3);
                /* p(2^j x) has roots cos(...) / 2^j */
                for (len = 1; len < p->length; len++)
                    fmpz_mul_2exp(p->coeffs + len, p->coeffs + len, j * len);
                nr = m;
                /* roots in (0, 1): all positive roots */
                *n01 = m / 2;
                break;

            case 2: /* sparse: x^m - 2 (a x - 1)^2 (close roots near 1/a) */
                /* keep this small since plain VCA is very slow here */
                m = 5 + n_randint(state, 40);
                fmpz_randtest_unsigned(a, state, 1 + n_randint(state, 8));
                fmpz_add_ui(a, a, 3);
                fmpz_poly_zero(p);
                fmpz_poly_set_coeff_si(p, m, 1);
                fmpz_poly_set_coeff_si(lin, 0, -1);
                fmpz_poly_set_coeff_fmpz(lin, 1, a);
                fmpz_poly_mul(q, lin, lin);
                fmpz_poly_scalar_mul_si(q, q, -2);
                fmpz_poly_add(p, p, q);
                /* not a known count in general; mark as unknown */
                nr = -1;
                break;

            case 3: /* root at +/- 2^e times a small polynomial: the root
                       can coincide with the power-of-two root bound */
                m = (slong) n_randint(state, 20) - 10;
                fmpz_poly_zero(lin);
                fmpz_poly_set_coeff_si(lin, 1, 1);
                fmpz_poly_set_coeff_si(lin, 0, n_randint(state, 2) ? 1 : -1);
                if (m >= 0)
                    fmpz_mul_2exp(lin->coeffs, lin->coeffs, m);
                else
                    fmpz_mul_2exp(lin->coeffs + 1, lin->coeffs + 1, -m);
                fmpz_poly_randtest_not_zero(p, state, 1 + n_randint(state, 6), 1 + n_randint(state, 4));
                fmpz_poly_mul(p, p, lin);
                nr = -1;
                break;

            default: /* generic random polynomial; count unknown */
                fmpz_poly_randtest_not_zero(p, state, 1 + n_randint(state, 60), 1 + n_randint(state, 200));
                nr = -1;
                break;
        }
    }
    while (fmpz_poly_is_zero(p) || !fmpz_poly_is_squarefree(p));

    fmpz_poly_clear(q);
    fmpz_poly_clear(lin);
    fmpz_clear(a);
    fmpz_clear(b);

    return nr;
}

TEST_FUNCTION_START(fmpz_poly_num_real_roots, state)
{
    int iter;

    /* call with random nonzero polynomials */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        slong k;
        fmpz_poly_t p;

        fmpz_poly_init(p);
        fmpz_poly_randtest_not_zero(p, state, 20, 10 + n_randint(state, 100));
        k = fmpz_poly_num_real_roots(p);
        if (k < 0 || k > fmpz_poly_degree(p))
        {
            printf("ERROR:\n");
            flint_printf("got k in wrong range k = %wd\n", k);
            printf("p = "); fmpz_poly_print(p); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(p);
    }

    for (iter = 0; iter < 5000 * flint_test_multiplier(); iter++)
    {
        slong k1, k2;
        fmpz_poly_t p;

        fmpz_poly_init(p);
        /* currently the code of num_real_roots only has a special branch */
        /* for length <= 5. We only test these cases.                     */
        do
        {
            fmpz_poly_randtest_not_zero(p, state, 1 + n_randint(state, 5), 100);
        } while (!fmpz_poly_is_squarefree(p));

        k1 = fmpz_poly_num_real_roots_sturm(p);
        k2 = fmpz_poly_num_real_roots(p);
        if (k1 != k2)
        {
            printf("ERROR:\n");
            flint_printf("found k1=%wd and k2=%wd\n", k1, k2);
            printf("p = "); fmpz_poly_print(p); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(p);
    }

    /* regression test: the root 2 coincides with the root bound */
    {
        fmpz_poly_t p;
        fmpz_poly_init(p);
        fmpz_poly_set_str(p, "4  -2 -1 -1 1");  /* (x - 2)(x^2 + x + 1) */
        if (fmpz_poly_num_real_roots_vca(p) != 1 ||
            fmpz_poly_num_real_roots_upper_bound(p) < 1)
        {
            flint_printf("FAIL: root at the root bound\n");
            flint_abort();
        }
        fmpz_poly_clear(p);
    }

    /* test the different algorithms against each other and against known
       root counts */
    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t p;
        slong nr, n01, k1, k2, k3, k4, k5, k6;
        int small, use_vca;

        fmpz_poly_init(p);
        nr = _randtest_known_roots(p, &n01, state);

        /* VCA (slow for some of these inputs, and tested separately) is
           only needed as a reference when the number of roots is not known */
        use_vca = (nr == -1);

        k1 = fmpz_poly_num_real_roots(p);
        k2 = use_vca ? fmpz_poly_num_real_roots_vca(p) : nr;
        /* avoid slow cases for the (subresultant) Sturm sequence */
        small = (p->length * FLINT_ABS(fmpz_poly_max_bits(p)) <= 3000);
        k3 = small ? fmpz_poly_num_real_roots_sturm(p) : k2;

        k4 = fmpz_poly_num_real_roots_0_1(p);
        k5 = use_vca ? fmpz_poly_num_real_roots_0_1_vca(p) : n01;
        k6 = small ? fmpz_poly_num_real_roots_0_1_sturm(p) : k5;

        if (k1 != k2 || k1 != k3 || (nr != -1 && k1 != nr) ||
            fmpz_poly_num_real_roots_upper_bound(p) < k1 ||
            k4 != k5 || k4 != k6 || (nr != -1 && k4 != n01))
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "); fmpz_poly_print(p); flint_printf("\n");
            flint_printf("expected %wd, %wd\n", nr, n01);
            flint_printf("num_real_roots: default %wd, vca %wd, sturm %wd\n", k1, k2, k3);
            flint_printf("num_real_roots_0_1: default %wd, vca %wd, sturm %wd\n", k4, k5, k6);
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(p);
    }

    TEST_FUNCTION_END(state);
}
