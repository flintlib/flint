/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "ulong_extras.h"

/*
    Targeted test of the deflation branch of fmpz_poly_factor and the
    Capelli-certification pathways inside it, comparing against
    fmpz_poly_factor_zassenhaus, which never takes the deflation branch
    and therefore serves as an independent oracle.

    The generators below are engineered to hit, in particular:
    - certified irreducible inflations, both via the exact rational p-th
      power test (linear deflation factors) and via witness search
      (p = 2 and odd p, including multiquadratic-type inputs whose
      witnesses have density 1/(p deg T));
    - reducible inflations with no witness (perfect-power theta, norm-form
      products g(x) g(-x)), exercising budget/evidence bails and fallback;
    - self-reciprocal (cyclotomic-guard) inputs;
    - chained composite deflations d = 4, 6, 8, 9, 12;
    - powers of x, contents, signs and multiplicities on top of deflation.
*/

static int
factor_matches(const fmpz_poly_factor_t a, const fmpz_poly_factor_t b)
{
    slong i, j;
    int * used;
    int ok = 1;

    if (a->num != b->num || !fmpz_equal(&a->c, &b->c))
        return 0;

    used = flint_calloc(b->num, sizeof(int));
    for (i = 0; i < a->num && ok; i++)
    {
        for (j = 0; j < b->num; j++)
        {
            if (!used[j] && a->exp[i] == b->exp[j]
                         && fmpz_poly_equal(a->p + i, b->p + j))
            {
                used[j] = 1;
                break;
            }
        }
        if (j == b->num)
            ok = 0;
    }
    flint_free(used);
    return ok;
}

static void
check_product(const fmpz_poly_factor_t fac, const fmpz_poly_t G,
              const char * which)
{
    slong i, j;
    fmpz_poly_t prod;

    fmpz_poly_init(prod);
    fmpz_poly_set_fmpz(prod, &fac->c);
    for (i = 0; i < fac->num; i++)
        for (j = 0; j < fac->exp[i]; j++)
            fmpz_poly_mul(prod, prod, fac->p + i);

    if (!fmpz_poly_equal(prod, G))
    {
        flint_printf("FAIL (%s: product)\n", which);
        flint_printf("G = "); fmpz_poly_print(G); flint_printf("\n");
        fflush(stdout);
        flint_abort();
    }
    fmpz_poly_clear(prod);
}

TEST_FUNCTION_START(fmpz_poly_factor_deflation_capelli, state)
{
    slong iter;

    for (iter = 0; iter < 150 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t G, T, t, u;
        fmpz_poly_factor_t fac, zfac;
        fmpz_t a;
        slong mode = n_randint(state, 6);
        ulong d;

        fmpz_poly_init(G); fmpz_poly_init(T);
        fmpz_poly_init(t); fmpz_poly_init(u);
        fmpz_init(a);

        switch (mode)
        {
        case 0:
            /* x^d - a and x^d + a, a sometimes a perfect power:
               exercises the exact rational p-th power test and odd-p
               witness search through chained deflation */
            {
                ulong ds[] = {4, 6, 8, 9, 12};
                d = ds[n_randint(state, 5)];
                fmpz_randtest_not_zero(a, state, 24);
                if (n_randint(state, 2))
                {
                    slong e = 2 + n_randint(state, 3);
                    fmpz_pow_ui(a, a, (ulong) e);   /* perfect power */
                }
                fmpz_poly_set_coeff_ui(G, (slong) d, 1);
                fmpz_neg(a, a);
                fmpz_poly_set_coeff_fmpz(G, 0, a);
            }
            break;

        case 1:
            /* +- g(x) g(-x): even, inflation of an irreducible piece is
               reducible with theta a square in K; no witness exists,
               exercising the evidence bail and fallback */
            do {
                fmpz_poly_randtest(t, state, 2 + n_randint(state, 4), 30);
            } while (fmpz_poly_degree(t) < 1 ||
                     fmpz_is_zero(t->coeffs + 0));
            {
                slong i;
                fmpz_poly_set(u, t);
                for (i = 1; i <= fmpz_poly_degree(u); i += 2)
                {
                    fmpz_poly_get_coeff_fmpz(a, u, i);
                    fmpz_neg(a, a);
                    fmpz_poly_set_coeff_fmpz(u, i, a);
                }
                fmpz_poly_mul(G, t, u);
            }
            break;

        case 2:
            /* minimal polynomial of sqrt(b) + sqrt(c):
               x^4 - 2(b+c) x^2 + (b-c)^2, multiquadratic-type input with
               rare witnesses; sometimes inflated once more */
            {
                slong b = 2 + n_randint(state, 30);
                slong c = 2 + n_randint(state, 30);
                fmpz_poly_set_coeff_ui(G, 4, 1);
                fmpz_set_si(a, -2 * (b + c));
                fmpz_poly_set_coeff_fmpz(G, 2, a);
                fmpz_set_si(a, (b - c) * (b - c));
                fmpz_poly_set_coeff_fmpz(G, 0, a);
                if (n_randint(state, 2))
                {
                    fmpz_poly_inflate(G, G, 2 + n_randint(state, 2));
                    if (fmpz_poly_degree(G) > 16)
                        fmpz_poly_deflate(G, G, 2);
                }
            }
            break;

        case 3:
            /* self-reciprocal, constant term 1: cyclotomic-guard path,
               e.g. x^(4k) + x^(2k) + 1 and x^(2k) + 1 */
            {
                slong k = 1 + n_randint(state, 4);
                if (n_randint(state, 2))
                {
                    fmpz_poly_set_coeff_ui(G, 4 * k, 1);
                    fmpz_poly_set_coeff_ui(G, 2 * k, 1);
                    fmpz_poly_set_coeff_ui(G, 0, 1);
                }
                else
                {
                    fmpz_poly_set_coeff_ui(G, 2 * k, 1);
                    fmpz_poly_set_coeff_ui(G, 0, 1);
                }
            }
            break;

        case 4:
            /* random polynomial in x^2 or x^3, occasionally sparse,
               coefficients up to ~60 bits: generic deflation exercise */
            do {
                fmpz_poly_randtest(t, state, 3 + n_randint(state, 6), 60);
            } while (fmpz_poly_degree(t) < 1);
            fmpz_poly_inflate(G, t, 2 + n_randint(state, 2));
            break;

        case 5:
            /* small multiquadratic norm form times an x^d - a piece:
               mixes certified and fallback levels in one chain */
            fmpz_poly_set_coeff_ui(t, 4, 1);
            fmpz_poly_set_coeff_si(t, 2, -10);
            fmpz_poly_set_coeff_ui(t, 0, 1);       /* min poly sqrt2+sqrt3 */
            fmpz_poly_set_coeff_ui(u, 4, 1);
            fmpz_randtest_not_zero(a, state, 16);
            fmpz_poly_set_coeff_fmpz(u, 0, a);
            fmpz_poly_mul(G, t, u);
            break;
        }

        if (fmpz_poly_degree(G) < 1)
            goto cleanup;

        /* dress the input: powers of x, content, sign, multiplicity */
        if (n_randint(state, 3) == 0)
            fmpz_poly_shift_left(G, G, 1 + n_randint(state, 3));
        if (n_randint(state, 3) == 0)
        {
            fmpz_randtest_not_zero(a, state, 10);
            fmpz_poly_scalar_mul_fmpz(G, G, a);
        }
        if (n_randint(state, 4) == 0 && fmpz_poly_degree(G) <= 12)
            fmpz_poly_mul(G, G, G);

        if (fmpz_poly_degree(G) > 40)
            goto cleanup;

        fmpz_poly_factor_init(fac);
        fmpz_poly_factor_init(zfac);

        fmpz_poly_factor(fac, G);
        fmpz_poly_factor_zassenhaus(zfac, G);

        check_product(fac, G, "factor");
        check_product(zfac, G, "zassenhaus");

        if (!factor_matches(fac, zfac))
        {
            flint_printf("FAIL (mismatch vs zassenhaus, mode %wd)\n", mode);
            flint_printf("G = "); fmpz_poly_print(G); flint_printf("\n");
            flint_printf("factor num = %wd, zassenhaus num = %wd\n",
                         fac->num, zfac->num);
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_factor_clear(fac);
        fmpz_poly_factor_clear(zfac);

cleanup:
        fmpz_poly_clear(G); fmpz_poly_clear(T);
        fmpz_poly_clear(t); fmpz_poly_clear(u);
        fmpz_clear(a);
    }

    TEST_FUNCTION_END(state);
}
