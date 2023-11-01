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
#define test_poly test_poly_roots
void test_poly(
    nmod_poly_factor_t roots,
    const nmod_poly_t f,
    int want_mult)
{
    slong i, multiplicity;
    nmod_poly_t q, qt, r;

    nmod_poly_init_mod(q, f->mod);
    nmod_poly_init_mod(qt, f->mod);
    nmod_poly_init_mod(r, f->mod);
    nmod_poly_set(q, f);

    nmod_poly_roots(roots, f, want_mult);

    for (i = 0; i < roots->num; i++)
    {
        if (nmod_poly_degree(roots->p + i) != 1)
        {
            flint_printf("FAIL:\ncheck root is linear\n");
            fflush(stdout);
            flint_abort();
        }

        if (roots->p[i].coeffs[1] != 1)
        {
            flint_printf("FAIL:\ncheck root is monic\n");
            fflush(stdout);
            flint_abort();
        }

        multiplicity = 0;
        while (nmod_poly_divrem(qt, r, q, roots->p + i), nmod_poly_is_zero(r))
        {
            nmod_poly_swap(q, qt);
            multiplicity++;
        }

        if (multiplicity <= 0)
        {
            flint_printf("FAIL:\ncheck root is a root\n");
            fflush(stdout);
            flint_abort();
        }

        if (roots->exp[i] != (want_mult ? multiplicity : 1))
        {
            flint_printf("FAIL:\ncheck root multiplicity\n");
            fflush(stdout);
            flint_abort();
        }
    }

    nmod_poly_roots(roots, q, want_mult);
    if (roots->num > 0)
    {
        flint_printf("FAIL:\ncheck missing roots\n");
        fflush(stdout);
        flint_abort();
    }

    nmod_poly_clear(q);
    nmod_poly_clear(qt);
    nmod_poly_clear(r);
}

TEST_FUNCTION_START(nmod_poly_factor_roots, state)
{
    slong i, j, k, l;

    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        mp_limb_t p;
        nmod_poly_t f;
        nmod_poly_factor_t r;

        p = n_randtest_prime(state, 1);

        nmod_poly_init(f, p);
        nmod_poly_factor_init(r);

        for (j = 0; j < 4; j++)
        {
            do {
                nmod_poly_randtest(f, state, n_randint(state, 25) + 1);
            } while (nmod_poly_is_zero(f));

            for (k = 0; k < 8; k++)
            {
                nmod_poly_t ff;
                nmod_poly_init(ff, p);
                nmod_poly_set_coeff_ui(ff, 1, 1);
                nmod_poly_set_coeff_ui(ff, 0, n_randint(state, p));
                for (l = 1 + n_randint(state, 8); l > 0; l--)
                    nmod_poly_mul(f, f, ff);
                nmod_poly_clear(ff);
            }

            if (n_randint(state, 2))
            {
                test_poly(r, f, 1);
                test_poly(r, f, 0);
            }
            else
            {
                test_poly(r, f, 0);
                test_poly(r, f, 1);
            }
        }

        nmod_poly_factor_clear(r);
        nmod_poly_clear(f);
    }

    TEST_FUNCTION_END(state);
}
#undef test_poly
