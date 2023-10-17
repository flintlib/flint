/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "nmod_poly_factor.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"

TEST_FUNCTION_START(fmpz_poly_hensel_start_continue_lift, state)
{
    int i, result;

    /* We check that lifting local factors of F yields factors */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t F, G, H, R;
        nmod_poly_factor_t f_fac;
        fmpz_poly_factor_t F_fac;
        slong bits, nbits, n, exp, j, part_exp;

        slong r;
        fmpz_poly_t *v, *w;
        slong *link;
        slong prev_exp;

        bits = n_randint(state, 200) + 1;
        nbits = n_randint(state, FLINT_BITS - 6) + 6;

        fmpz_poly_init(F);
        fmpz_poly_init(G);
        fmpz_poly_init(H);
        fmpz_poly_init(R);
        nmod_poly_factor_init(f_fac);
        fmpz_poly_factor_init(F_fac);

        n = n_randprime(state, nbits, 0);
        exp = bits / (FLINT_BIT_COUNT(n) - 1) + 1;
        part_exp = n_randint(state, exp);

        /* Produce F as the product of random G and H */
        {
            nmod_poly_t f;

            nmod_poly_init(f, n);

            do {
                do {
                    fmpz_poly_randtest(G, state, n_randint(state, 200) + 2, bits);
                } while (G->length < 2);

                fmpz_randtest_not_zero(G->coeffs, state, bits);
                fmpz_one(fmpz_poly_lead(G));

                do {
                    fmpz_poly_randtest(H, state, n_randint(state, 200) + 2, bits);
                } while (H->length < 2);

                fmpz_randtest_not_zero(H->coeffs, state, bits);
                fmpz_one(fmpz_poly_lead(H));

                fmpz_poly_mul(F, G, H);

                fmpz_poly_get_nmod_poly(f, F);
            } while (!nmod_poly_is_squarefree(f));

            fmpz_poly_get_nmod_poly(f, G);
            nmod_poly_factor_insert(f_fac, f, 1);
            fmpz_poly_get_nmod_poly(f, H);
            nmod_poly_factor_insert(f_fac, f, 1);
            nmod_poly_clear(f);
        }

        r = f_fac->num;
        v = flint_malloc((2*r - 2)*sizeof(fmpz_poly_t));
        w = flint_malloc((2*r - 2)*sizeof(fmpz_poly_t));
        link = flint_malloc((2*r - 2)*sizeof(slong));

        for (j = 0; j < 2*r - 2; j++)
        {
            fmpz_poly_init(v[j]);
            fmpz_poly_init(w[j]);
        }

        if (part_exp < 1)
        {
            _fmpz_poly_hensel_start_lift(F_fac, link, v, w, F, f_fac, exp);
        }
        else
        {
            fmpz_t nn;

            fmpz_init_set_ui(nn, n);

            prev_exp = _fmpz_poly_hensel_start_lift(F_fac, link, v, w,
                F, f_fac, part_exp);
            _fmpz_poly_hensel_continue_lift(F_fac, link, v, w,
                F, prev_exp, part_exp, exp, nn);

            fmpz_clear(nn);
        }

        result = 1;
        for (j = 0; j < F_fac->num; j++)
        {
            fmpz_poly_rem(R, F, F_fac->p + j);
            result &= (R->length == 0);
        }

        for (j = 0; j < 2*r - 2; j++)
        {
            fmpz_poly_clear(v[j]);
            fmpz_poly_clear(w[j]);
        }

        flint_free(link);
        flint_free(v);
        flint_free(w);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("bits = %wd, n = %wd, exp = %wd\n", bits, n, exp);
            fmpz_poly_print(F); flint_printf("\n\n");
            fmpz_poly_print(G); flint_printf("\n\n");
            fmpz_poly_print(H); flint_printf("\n\n");
            fmpz_poly_factor_print(F_fac); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_factor_clear(f_fac);
        fmpz_poly_factor_clear(F_fac);

        fmpz_poly_clear(F);
        fmpz_poly_clear(H);
        fmpz_poly_clear(G);
        fmpz_poly_clear(R);
    }

    TEST_FUNCTION_END(state);
}
