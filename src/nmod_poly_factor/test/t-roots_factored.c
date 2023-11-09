/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_poly.h"
#include "nmod_poly_factor.h"

/* Defined in t-roots.c and t-roots_factored.c */
#define test_poly test_poly_factored
void test_poly(
    nmod_poly_factor_t roots,
    const nmod_poly_t f,
    int want_mult,
    const n_factor_t * n)
{
    slong i, multiplicity;
    nmod_poly_t q, qt, r;

    nmod_poly_init_mod(q, f->mod);
    nmod_poly_init_mod(qt, f->mod);
    nmod_poly_init_mod(r, f->mod);

    if (!nmod_poly_roots_factored(roots, f, want_mult, n))
    {
        flint_printf("FAILED:\ncheck roots could be computed\n");
        fflush(stdout);
        flint_abort();
    }

    for (i = 0; i < roots->num; i++)
    {
        if (nmod_poly_degree(roots->p + i) != 1)
        {
            flint_printf("FAILED:\ncheck root is linear\n");
            fflush(stdout);
            flint_abort();
        }

        if (roots->p[i].coeffs[1] != 1)
        {
            flint_printf("FAILED:\ncheck root is monic\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_set(q, f);

        multiplicity = 0;
        while (nmod_poly_divrem(qt, r, q, roots->p + i), nmod_poly_is_zero(r))
        {
            nmod_poly_swap(q, qt);
            multiplicity++;
        }

        if (multiplicity <= 0)
        {
            flint_printf("FAILED:\ncheck root is a root\n");
            fflush(stdout);
            flint_abort();
        }

        if (roots->exp[i] != (want_mult ? multiplicity : 1))
        {
            flint_printf("FAILED:\ncheck root multiplicity\n");
            fflush(stdout);
            flint_abort();
        }
    }

    if (f->mod.n < 4000)
    {
        ulong k;

        for (k = 0; k < f->mod.n; k++)
        {
            int found = 0;

            if (0 != nmod_poly_evaluate_nmod(f, k))
                continue;

            for (i = 0; i < roots->num; i++)
            {
                if (0 == nmod_poly_evaluate_nmod(roots->p + i, k))
                {
                    if (found)
                    {
                        flint_printf("FAILED:\ncheck duplicate roots\n");
                        fflush(stdout);
                        flint_abort();
                    }
                    found = 1;
                }
            }

            if (!found)
            {
                flint_printf("FAILED:\ncheck missing roots\n");
                fflush(stdout);
                flint_abort();
            }
        }
    }

    nmod_poly_clear(q);
    nmod_poly_clear(qt);
    nmod_poly_clear(r);
}

TEST_FUNCTION_START(nmod_poly_factor_roots_factored, state)
{
    slong i, j, k, l;

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f;
        nmod_poly_factor_t roots;
        mp_limb_t a, n;
        mp_limb_t * sqrt;
        n_factor_t nfac;

        n = n_randtest_bits(state, n_randint(state, FLINT_BITS) + 1);
        n = FLINT_MAX(n, UWORD(2));
        n_factor_init(&nfac);
        n_factor(&nfac, n, 0);

        nmod_poly_init(f, n);
        nmod_poly_factor_init(roots);

        for (j = 0; j < 50; j++)
        {
            a = n_randint(state, n);

            nmod_poly_zero(f);
            nmod_poly_set_coeff_ui(f, 2, n - 1);
            nmod_poly_set_coeff_ui(f, 0, a);

            if (!nmod_poly_roots_factored(roots, f, 0, &nfac))
            {
                flint_printf("FAILED:\ncheck sqrt could be calculated\n");
                fflush(stdout);
                flint_abort();
            }

            if (roots->num != n_sqrtmodn(&sqrt, a, &nfac))
            {
                flint_printf("FAILED:\ncheck root count against n_sqrtmodn\n");
                fflush(stdout);
                flint_abort();
            }

            if (sqrt)
                flint_free(sqrt);

            test_poly(roots, f, 0, &nfac);
            test_poly(roots, f, 1, &nfac);
        }

        nmod_poly_clear(f);
        nmod_poly_factor_clear(roots);
    }

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        ulong p;
        nmod_poly_t f;
        nmod_poly_factor_t r;
        n_factor_t n;

        p = n_randbits(state, n_randint(state, FLINT_BITS) + 1);
        p = FLINT_MAX(p, UWORD(2));

        n_factor_init(&n);
        n_factor(&n, p, 1);

        nmod_poly_init(f, p);
        nmod_poly_factor_init(r);

        for (j = 0; j < 4; j++)
        {
            slong m = 80/FLINT_BIT_COUNT(p);

            do {
                nmod_poly_randtest(f, state, n_randint(state, 10 + m) + 1);
            } while (nmod_poly_is_zero(f));

            for (k = 0; k < m; k++)
            {
                nmod_poly_t ff;
                nmod_poly_init_mod(ff, f->mod);
                nmod_poly_set_coeff_ui(ff, 1, 1);
                nmod_poly_set_coeff_ui(ff, 0, n_randint(state, p));
                for (l = n_randint(state, m); l > 0; l--)
                    nmod_poly_mul(f, f, ff);
                nmod_poly_clear(ff);
            }

            test_poly(r, f, 0, &n);
            if (r->num < 1000)
                test_poly(r, f, 1, &n);
        }

        nmod_poly_factor_clear(r);
        nmod_poly_clear(f);
    }

    TEST_FUNCTION_END(state);
}
#undef test_poly
