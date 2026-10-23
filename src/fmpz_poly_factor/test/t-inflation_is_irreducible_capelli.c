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

/* is G irreducible, decided by zassenhaus (never uses deflation) */
static int
zassenhaus_irreducible(const fmpz_poly_t G)
{
    fmpz_poly_factor_t fac;
    int res;
    fmpz_poly_factor_init(fac);
    fmpz_poly_factor_zassenhaus(fac, G);
    res = (fac->num == 1 && fac->exp[0] == 1
                         && fmpz_poly_degree(fac->p + 0)
                            == fmpz_poly_degree(G));
    fmpz_poly_factor_clear(fac);
    return res;
}

TEST_FUNCTION_START(fmpz_poly_factor_inflation_is_irreducible_capelli, state)
{
    slong iter;

    /* trivial and edge inputs */
    {
        fmpz_poly_t T;
        fmpz_poly_init(T);

        fmpz_poly_set_ui(T, 5);                       /* deg < 1 */
        FLINT_TEST(!_fmpz_poly_factor_inflation_is_irreducible_capelli(T, 2));
        FLINT_TEST(!fmpz_poly_factor_inflation_is_irreducible_capelli(T, 2));

        fmpz_poly_set_coeff_ui(T, 1, 1);
        fmpz_poly_set_coeff_ui(T, 0, 0);              /* T = x: theta = 0 */
        FLINT_TEST(!_fmpz_poly_factor_inflation_is_irreducible_capelli(T, 2));

        fmpz_poly_set_coeff_si(T, 0, 1);              /* T = x + 1, d = 0, 1 */
        FLINT_TEST(!fmpz_poly_factor_inflation_is_irreducible_capelli(T, 0));
        FLINT_TEST(fmpz_poly_factor_inflation_is_irreducible_capelli(T, 1));

        fmpz_poly_clear(T);
    }

    /* degree 1: the exact rational p-th power test, all branches,
       including non-primitive input and negative leading coefficient */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t T, G;
        fmpz_t a, b, g;
        ulong ps[] = {2, 3, 5};
        ulong p = ps[n_randint(state, 3)];
        int certified;

        fmpz_poly_init(T); fmpz_poly_init(G);
        fmpz_init(a); fmpz_init(b); fmpz_init(g);

        fmpz_randtest_not_zero(a, state, 40);
        if (n_randint(state, 2))
            fmpz_pow_ui(a, a, p);                     /* perfect power */
        do {
            fmpz_randtest_not_zero(b, state, 8);
        } while (fmpz_is_zero(b));
        if (n_randint(state, 2))
            fmpz_one(b);
        if (n_randint(state, 3) == 0)                 /* force common factor */
        {
            fmpz_mul(a, a, b);
        }

        /* T = b x - a, theta = a / b; sometimes negate the whole thing */
        fmpz_poly_set_coeff_fmpz(T, 1, b);
        fmpz_neg(a, a);
        fmpz_poly_set_coeff_fmpz(T, 0, a);
        if (n_randint(state, 2))
            fmpz_poly_neg(T, T);

        /* skip if T is not irreducible, i.e. degenerate */
        if (fmpz_poly_degree(T) != 1)
            goto cleanup1;

        /* the public function normalises content and sign; also exercise
           the underscore version on the normalised polynomial */
        certified = fmpz_poly_factor_inflation_is_irreducible_capelli(T, p);
        {
            fmpz_t cont;
            fmpz_init(cont);
            fmpz_poly_content(cont, T);
            if (fmpz_sgn(fmpz_poly_lead(T)) < 0)
                fmpz_neg(cont, cont);
            fmpz_poly_scalar_divexact_fmpz(G, T, cont);
            fmpz_clear(cont);
            if (certified !=
                _fmpz_poly_factor_inflation_is_irreducible_capelli(G, p))
            {
                flint_printf("FAIL (deg 1 normalisation, p = %wu)\n", p);
                flint_abort();
            }
        }
        fmpz_poly_inflate(G, T, p);

        /* the linear case is decided exactly: certified iff irreducible */
        if (certified != zassenhaus_irreducible(G))
        {
            flint_printf("FAIL (deg 1, p = %wu)\nT = ", p);
            fmpz_poly_print(T); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

cleanup1:
        fmpz_poly_clear(T); fmpz_poly_clear(G);
        fmpz_clear(a); fmpz_clear(b); fmpz_clear(g);
    }

    /* degree 2 with p = 2: the exact quadratic-field square test.
       Generate both random irreducible quadratics and engineered
       norm-form reducible inflations a(x^2+bx+c)(x^2-bx+c). */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t T, G, u, v;
        int certified;

        fmpz_poly_init(T); fmpz_poly_init(G);
        fmpz_poly_init(u); fmpz_poly_init(v);

        if (n_randint(state, 2))
        {
            fmpz_poly_randtest(T, state, 3, 12);
        }
        else
        {
            /* T = deflation of (x^2+bx+c)(x^2-bx+c): theta is a square */
            slong b = 1 + (slong) n_randint(state, 20);
            slong c = (slong) n_randint(state, 40) - 20;
            fmpz_poly_set_coeff_ui(u, 2, 1);
            fmpz_poly_set_coeff_si(u, 1, b);
            fmpz_poly_set_coeff_si(u, 0, c);
            fmpz_poly_set_coeff_ui(v, 2, 1);
            fmpz_poly_set_coeff_si(v, 1, -b);
            fmpz_poly_set_coeff_si(v, 0, c);
            fmpz_poly_mul(G, u, v);
            fmpz_poly_deflate(T, G, 2);
        }

        if (fmpz_poly_degree(T) != 2 || !zassenhaus_irreducible(T))
            goto cleanup2;

        certified = _fmpz_poly_factor_inflation_is_irreducible_capelli(T, 2);
        fmpz_poly_inflate(G, T, 2);

        /* the quadratic case with p = 2 is decided exactly */
        if (certified != zassenhaus_irreducible(G))
        {
            flint_printf("FAIL (deg 2, p = 2)\nT = ");
            fmpz_poly_print(T); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

cleanup2:
        fmpz_poly_clear(T); fmpz_poly_clear(G);
        fmpz_poly_clear(u); fmpz_poly_clear(v);
    }

    /* general degrees: soundness (certified implies irreducible) for the
       witness search, over irreducible factors of random and structured
       polynomials, prime and composite d */
    for (iter = 0; iter < 40 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t t, T, G;
        fmpz_poly_factor_t fac;
        ulong ds[] = {2, 3, 4, 5, 6, 9};
        ulong d = ds[n_randint(state, 6)];
        int certified;

        fmpz_poly_init(t); fmpz_poly_init(T); fmpz_poly_init(G);
        fmpz_poly_factor_init(fac);

        switch (n_randint(state, 3))
        {
        case 0:                                   /* random */
            fmpz_poly_randtest(t, state, 3 + n_randint(state, 5), 20);
            break;
        case 1:                                   /* multiquadratic-type */
            {
                slong b = 2 + (slong) n_randint(state, 20);
                slong c = 2 + (slong) n_randint(state, 20);
                fmpz_poly_set_coeff_ui(t, 4, 1);
                fmpz_poly_set_coeff_si(t, 2, -2 * (b + c));
                fmpz_poly_set_coeff_si(t, 0, (b - c) * (b - c));
            }
            break;
        case 2:                                   /* self-reciprocal */
            fmpz_poly_set_coeff_ui(t, 4, 1);
            fmpz_poly_set_coeff_si(t, 3, -1);
            fmpz_poly_set_coeff_ui(t, 2, 1);
            fmpz_poly_set_coeff_si(t, 1, -1);
            fmpz_poly_set_coeff_ui(t, 0, 1);
            break;
        }

        if (fmpz_poly_degree(t) < 1)
            goto cleanup3;

        fmpz_poly_factor(fac, t);
        if (fac->num < 1 || fmpz_poly_degree(fac->p + 0) < 1)
            goto cleanup3;
        fmpz_poly_set(T, fac->p + 0);

        if ((ulong) fmpz_poly_degree(T) * d > 24)
            goto cleanup3;                        /* keep the zassenhaus cheap */

        certified = fmpz_poly_factor_inflation_is_irreducible_capelli(T, d);
        fmpz_poly_inflate(G, T, d);

        if (certified && !zassenhaus_irreducible(G))
        {
            flint_printf("FAIL (soundness, d = %wu)\nT = ", d);
            fmpz_poly_print(T); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

cleanup3:
        fmpz_poly_factor_clear(fac);
        fmpz_poly_clear(t); fmpz_poly_clear(T); fmpz_poly_clear(G);
    }

    /* deterministic large abelian case: the degree-64 deflation T of the
       Swinnerton-Dyer polynomial S_7 (minimal polynomial of
       sqrt2+...+sqrt17).  This reaches the phase-B specialist search:
       witnesses exist only at the density-2^-7 primes where all seven
       p_i are quadratic nonresidues, in degree-1 residue fields.  The
       prime scan is deterministic, so the certificate is reproducible. */
    {
        fmpz_poly_t g, T;

        fmpz_poly_init(g);
        fmpz_poly_init(T);

        fmpz_poly_swinnerton_dyer(g, 7);
        fmpz_poly_deflate(T, g, 2);

        FLINT_TEST(_fmpz_poly_factor_inflation_is_irreducible_capelli(T, 2));

        fmpz_poly_clear(g);
        fmpz_poly_clear(T);
    }

    /* odd-p phase B: with p = 61 and T = y^2 - 2 the inflated degree
       exceeds the phase-B threshold, phase A carries no Euler information
       (61 rarely divides q^f - 1 there), and phase B decides via primes
       q = 1 mod 61.  sqrt(2) has norm -2, not a 61st power, so x^61 -
       sqrt(2) is irreducible over Q(sqrt 2); the scan is deterministic. */
    {
        fmpz_poly_t T;
        fmpz_poly_init(T);
        fmpz_poly_set_coeff_ui(T, 2, 1);
        fmpz_poly_set_coeff_si(T, 0, -2);
        FLINT_TEST(_fmpz_poly_factor_inflation_is_irreducible_capelli(T, 61));
        fmpz_poly_clear(T);
    }

    TEST_FUNCTION_END(state);
}
