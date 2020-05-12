/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mod_poly_factor.h"
#include "ulong_extras.h"


void test_poly(
    fmpz_mod_poly_factor_t roots,
    const fmpz_mod_poly_t f,
    int want_mult,
    const fmpz_factor_t n)
{
    slong i, multiplicity;
    fmpz_mod_poly_t q, qt, r;

    fmpz_mod_poly_init(q, &f->p);
    fmpz_mod_poly_init(qt, &f->p);
    fmpz_mod_poly_init(r, &f->p);

    if (!fmpz_mod_poly_roots_factored(roots, f, want_mult, n))
    {
        flint_printf("FAILED:\ncheck roots could be computed\n");
        flint_abort();
    }

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

        fmpz_mod_poly_set(q, f);

        multiplicity = 0;
        while (fmpz_mod_poly_divrem(qt, r, q, roots->poly + i), fmpz_mod_poly_is_zero(r))
        {
            fmpz_mod_poly_swap(q, qt);
            multiplicity++;
        }

        if (multiplicity <= 0)
        {
            flint_printf("FAILED:\ncheck root is a root\n");
            fmpz_mod_poly_print_pretty(roots->poly + i, "x"); printf("\n");
            flint_abort();
        }

        if (roots->exp[i] != (want_mult ? multiplicity : 1))
        {
            flint_printf("FAILED:\ncheck root multiplicity\n");
            flint_abort();
        }
    }

    if (fmpz_cmp_si(&f->p, 1000) < 0)
    {
        fmpz_t e, k;

        fmpz_init(e);
        fmpz_init(k);

        for (fmpz_zero(k); fmpz_cmp(k, &f->p) < 0; fmpz_add_ui(k, k, 1))
        {
            int found = 0;

            fmpz_mod_poly_evaluate_fmpz(e, f, k);
            if (!fmpz_is_zero(e))
                continue;

            for (i = 0; i < roots->num; i++)
            {
                fmpz_mod_poly_evaluate_fmpz(e, roots->poly + i, k);
                if (fmpz_is_zero(e))
                {
                    if (found)
                    {
                        flint_printf("FAILED:\ncheck duplicate roots\n");
                        flint_abort();
                    }
                    found = 1;
                }
            }

            if (!found)
            {
                flint_printf("FAILED:\ncheck missing roots\n");
                flint_abort();
            }
        }

        fmpz_clear(k);
        fmpz_clear(e);
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

    flint_printf("roots_factored....");
    fflush(stdout);

    {
        fmpz_t p, p2;
        fmpz_mod_poly_t f;
        fmpz_mod_poly_factor_t r;
        fmpz_factor_t n;

        fmpz_init_set_ui(p, UWORD_MAX_PRIME);
        fmpz_init(p2);
        fmpz_pow_ui(p2, p, 2);
        fmpz_mod_poly_init(f, p2);
        fmpz_mod_poly_factor_init(r);
        fmpz_factor_init(n);
        _fmpz_factor_append(n, p, 2);

        fmpz_mod_poly_set_coeff_fmpz(f, 1, p);

        if (fmpz_mod_poly_roots_factored(r, f, 0, n))
        {
            flint_printf("FAILED:\ncheck non example with too many roots\n");
            flint_abort();
        }

        fmpz_factor_clear(n);

        fmpz_mod_poly_factor_clear(r);
        fmpz_mod_poly_clear(f);
        fmpz_clear(p);
        fmpz_clear(p2);
    }

    {
        fmpz_t one, p, q, p2q;
        fmpz_mod_poly_t f;
        fmpz_mod_poly_factor_t r;
        fmpz_factor_t n;
        ulong tp = n_nextprime(UWORD(1) << (FLINT_BITS - 1), 1);
        ulong tq = n_nextprime(tp, 1);

        fmpz_init_set_ui(one, 1);
        fmpz_init_set_ui(p, tp);
        fmpz_init_set_ui(q, tq);
        fmpz_init_set(p2q, q);
        fmpz_mul(p2q, p2q, p);
        fmpz_mul(p2q, p2q, p);

        fmpz_mod_poly_init(f, p2q);
        fmpz_mod_poly_factor_init(r);
        fmpz_factor_init(n);
        _fmpz_factor_append(n, p, 2);
        _fmpz_factor_append(n, q, 1);

        fmpz_mod_poly_set_coeff_fmpz(f, 0, one);
        fmpz_mod_poly_set_coeff_fmpz(f, 1, q);
        fmpz_mod_poly_scalar_mul_fmpz(f, f, p);

        if ((!fmpz_mod_poly_roots_factored(r, f, 0, n)) || (r->num != 0))
        {
            flint_printf("FAILED:\ncheck example with no roots\n");
            flint_abort();
        }

        fmpz_factor_clear(n);
        fmpz_mod_poly_factor_clear(r);
        fmpz_mod_poly_clear(f);
        fmpz_clear(p);
        fmpz_clear(q);
        fmpz_clear(p2q);
        fmpz_clear(one);
    }

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t f;
        fmpz_mod_poly_factor_t r;
        fmpz_factor_t n;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 9);
        fmpz_add_ui(p, p, 2);
        for (j = n_randint(state, 3) + 1; j > 0; j--)
        {
            fmpz_t q;
            fmpz_init(q);
            fmpz_randtest_unsigned(q, state, 6);
            fmpz_add_ui(q, q, 3);
            fmpz_nextprime(q, q, 1);
            fmpz_pow_ui(q, q, n_randint(state, 3) + 1);
            fmpz_mul(p, p, q);
        }

        fmpz_factor_init(n);
        fmpz_factor(n, p);

        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_factor_init(r);

        for (j = 0; j < 4; j++)
        {
            slong m = 80/fmpz_bits(p);

            do {
                fmpz_mod_poly_randtest(f, state, n_randint(state, 6 + m) + 1);
            } while (fmpz_mod_poly_is_zero(f));

            for (k = 0; k < m; k++)
            {
                fmpz_mod_poly_t ff;
                fmpz_mod_poly_init(ff, p);
                fmpz_mod_poly_fit_length(ff, 2);
                fmpz_one(ff->coeffs + 1);
                fmpz_randm(ff->coeffs + 0, state, p);
                ff->length = 2;
                for (l = n_randint(state, m); l > 0; l--)
                    fmpz_mod_poly_mul(f, f, ff);
                fmpz_mod_poly_clear(ff);
            }

            test_poly(r, f, 0, n);
            if (r->num < 1000)
                test_poly(r, f, 1, n);
        }

        fmpz_mod_poly_factor_clear(r);
        fmpz_mod_poly_clear(f);
        fmpz_factor_clear(n);
        fmpz_clear(p);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
