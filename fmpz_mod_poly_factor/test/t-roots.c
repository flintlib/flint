/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mod_poly_factor.h"


void test_poly(
    fmpz_mod_poly_factor_t roots,
    const fmpz_mod_poly_t f,
    int want_mult)
{
    slong i, multiplicity;
    fmpz_mod_poly_t q, qt, r;

    fmpz_mod_poly_init(q, &f->p);
    fmpz_mod_poly_init(qt, &f->p);
    fmpz_mod_poly_init(r, &f->p);
    fmpz_mod_poly_set(q, f);

    fmpz_mod_poly_roots(roots, f, want_mult);

    for (i = 0; i < roots->num; i++)
    {
        if (fmpz_mod_poly_degree(roots->poly + i) != 1)
        {
            flint_printf("FAILED:\ncheck root is linear\n");
            flint_abort();
        }

        if (!fmpz_is_one(roots->poly[i].coeffs + 1) != 0)
        {
            flint_printf("FAILED:\ncheck root is monic\n");
            flint_abort();
        }

        multiplicity = 0;
        while (fmpz_mod_poly_divrem(qt, r, q, roots->poly + i), fmpz_mod_poly_is_zero(r))
        {
            fmpz_mod_poly_swap(q, qt);
            multiplicity++;
        }

        if (multiplicity <= 0)
        {
            flint_printf("FAILED:\ncheck root is a root\n");
            flint_abort();
        }

        if (roots->exp[i] != (want_mult ? multiplicity : 1))
        {
            flint_printf("FAILED:\ncheck root multiplicity\n");
            flint_abort();
        }
    }

    fmpz_mod_poly_roots(roots, q, want_mult);
    if (roots->num > 0)
    {
        flint_printf("FAILED:\ncheck missing roots\n");
        flint_abort();
    }

    fmpz_mod_poly_clear(q);
    fmpz_mod_poly_clear(qt);
    fmpz_mod_poly_clear(r);
}


int
main(void)
{
    slong i, j, k, l;
    FLINT_TEST_INIT(state);

    flint_printf("roots....");
    fflush(stdout);

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t f;
        fmpz_mod_poly_factor_t r;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 100);

        fmpz_nextprime(p, p, 1);
        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_factor_init(r);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpz_mod_poly_randtest(f, state, n_randint(state, 20) + 1);
            } while (fmpz_mod_poly_is_zero(f));

            for (k = 0; k < 5; k++)
            {
                fmpz_mod_poly_t ff;
                fmpz_mod_poly_init(ff, p);
                fmpz_mod_poly_fit_length(ff, 2);
                fmpz_one(ff->coeffs + 1);
                fmpz_randm(ff->coeffs + 0, state, p);
                ff->length = 2;
                for (l = 1 + n_randint(state, 5); l > 0; l--)
                    fmpz_mod_poly_mul(f, f, ff);
                fmpz_mod_poly_clear(ff);
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

        fmpz_mod_poly_factor_clear(r);
        fmpz_mod_poly_clear(f);
        fmpz_clear(p);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
