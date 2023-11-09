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

TEST_FUNCTION_START(fmpz_poly_hensel_lift, state)
{
    int i, result;

    /* We check that lifting local factors of F_poly yields factors */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t F_poly, F_poly2, F_poly3, A, B, G, H,
            A_out, B_out, G_out, H_out, Prod_1, Prod_2;
        nmod_poly_t a, b, d, g, h, prod;
        fmpz_t p, p1, big_P, p1_2, big_P_2;
        slong bits, length, nbits, n, exp, part_exp;

        bits = n_randint(state, 200) + 1;
        nbits = n_randint(state, FLINT_BITS - 6) + 6;

        fmpz_init(p);
        fmpz_init(p1);
        fmpz_init(big_P);
        fmpz_init(p1_2);
        fmpz_init(big_P_2);

        fmpz_poly_init(F_poly);
        fmpz_poly_init(F_poly2);
        fmpz_poly_init(F_poly3);

        fmpz_poly_init(Prod_1);
        fmpz_poly_init(Prod_2);

        fmpz_poly_init(A);
        fmpz_poly_init(B);
        fmpz_poly_init(G);
        fmpz_poly_init(H);
        fmpz_poly_init(A_out);
        fmpz_poly_init(B_out);
        fmpz_poly_init(G_out);
        fmpz_poly_init(H_out);

        n = n_randprime(state, nbits, 0);
        exp = bits/(FLINT_BIT_COUNT(n) - 1) + 1;
        part_exp = n_randint(state, exp);

        nmod_poly_init(g, n);
        nmod_poly_init(h, n);
        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(d, n);
        nmod_poly_init(prod, n);

        do {
            length = n_randint(state, 200) + 2;

            do { fmpz_poly_randtest(F_poly2, state, length, bits); }
            while (F_poly2->length < 2);

            fmpz_set_ui(F_poly2->coeffs, n_randbits(state, FLINT_MIN(bits, SMALL_FMPZ_BITCOUNT_MAX)));
            fmpz_set_ui(F_poly2->coeffs + F_poly2->length - 1, 1);

            length = n_randint(state, 200) + 2;

            do { fmpz_poly_randtest(F_poly3, state, length, bits); }
            while (F_poly3->length < 2);

            fmpz_set_ui(F_poly3->coeffs, n_randbits(state, FLINT_MIN(bits, SMALL_FMPZ_BITCOUNT_MAX)));
            fmpz_set_ui(F_poly3->coeffs + F_poly3->length - 1, 1);

            fmpz_poly_mul(F_poly, F_poly2, F_poly3);

            fmpz_poly_get_nmod_poly(prod, F_poly);
        } while (!nmod_poly_is_squarefree(prod));

        fmpz_poly_get_nmod_poly(g, F_poly2);
        fmpz_poly_get_nmod_poly(h, F_poly3);

        nmod_poly_xgcd(d, a, b, g, h);
        nmod_poly_clear(prod);
        nmod_poly_clear(d);

        fmpz_poly_set_nmod_poly(A, a);
        fmpz_poly_set_nmod_poly(B, b);
        fmpz_poly_set_nmod_poly(G, g);
        fmpz_poly_set_nmod_poly(H, h);

        fmpz_set_ui(p, n);
        fmpz_set_ui(p1, n);
        fmpz_set_ui(big_P, n);
        fmpz_set_ui(p1_2, n);
        fmpz_set_ui(big_P_2, n);

        part_exp = 1;

        while (part_exp < exp)
        {
            fmpz_set(p, big_P);
            fmpz_set_ui(p1, n);
            fmpz_set_ui(big_P, n);

            if (exp - part_exp <= part_exp)
            {
                fmpz_pow_ui(p1, p1, exp - part_exp);
                fmpz_pow_ui(big_P, big_P, exp);
                part_exp = exp;
            }
            else
            {
                fmpz_set(p1, p);
                fmpz_pow_ui(big_P, big_P, 2*part_exp);
                part_exp = 2*part_exp;
            }

            fmpz_poly_hensel_lift(G_out, H_out, A_out, B_out, F_poly,
                G, H, A, B, p, p1);

            fmpz_poly_set(G, G_out);
            fmpz_poly_set(H, H_out);
            fmpz_poly_set(A, A_out);
            fmpz_poly_set(B, B_out);
        }

        fmpz_poly_mul(Prod_1, A, G);
        fmpz_poly_mul(Prod_2, B, H);
        fmpz_poly_add(Prod_1, Prod_1, Prod_2);

        fmpz_poly_scalar_smod_fmpz(Prod_1, Prod_1, big_P);

        result = (Prod_1->length == 1 && fmpz_is_one(Prod_1->coeffs));

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("length = %wd, bits = %wd, n = %wd, exp = %wd\n", length, bits, n, exp);
            fmpz_poly_print(F_poly); flint_printf("\n\n");
            fmpz_poly_print(F_poly2); flint_printf("\n\n");
            fmpz_poly_print(F_poly3); flint_printf("\n\n");
            fmpz_poly_print(Prod_1); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(g);
        nmod_poly_clear(h);
        nmod_poly_clear(a);
        nmod_poly_clear(b);

        fmpz_poly_clear(Prod_1);
        fmpz_poly_clear(Prod_2);

        fmpz_poly_clear(A);
        fmpz_poly_clear(B);
        fmpz_poly_clear(G);
        fmpz_poly_clear(H);
        fmpz_poly_clear(A_out);
        fmpz_poly_clear(B_out);
        fmpz_poly_clear(G_out);
        fmpz_poly_clear(H_out);

        fmpz_clear(p);
        fmpz_clear(p1);
        fmpz_clear(big_P);
        fmpz_clear(p1_2);
        fmpz_clear(big_P_2);

        fmpz_poly_clear(F_poly3);
        fmpz_poly_clear(F_poly2);
        fmpz_poly_clear(F_poly);
    }

    TEST_FUNCTION_END(state);
}
