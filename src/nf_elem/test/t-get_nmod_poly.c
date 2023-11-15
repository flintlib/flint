/*
    Copyright (C) 2018 Tommy Hofmann
                  2020 Julian Rüth

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_poly.h"
#include "nf_elem.h"

TEST_FUNCTION_START(nf_elem_get_nmod_poly, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong j;
        nf_t nf;
        nf_elem_t a;
        nmod_poly_t reduced_elem;
        ulong mod;
        fmpz_t coeff;

        mod = n_randtest_not_zero(state);

        fmpz_init(coeff);
        nmod_poly_init(reduced_elem, mod);

        nf_init_randtest(nf, state, 40, 200);

        nf_elem_init(a, nf);

        nf_elem_randtest(a, state, 200, nf);

        nf_elem_get_nmod_poly_den(reduced_elem, a, nf, 0);

        for (j = 0; j < fmpq_poly_degree(nf->pol); j++)
        {
            nf_elem_get_coeff_fmpz(coeff, a, j, nf);
            result = (nmod_poly_get_coeff_ui(reduced_elem, j) == fmpz_fdiv_ui(coeff, mod));
            if (!result)
            {
                printf("FAIL: Reducing without denominator\n");
                printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
                printf("n = "); flint_printf("%u\n", mod);
                printf("a mod n = "); nmod_poly_print_pretty(reduced_elem, "x"); printf("\n");
                flint_abort();
            }
        }

        nf_elem_clear(a, nf);
        nmod_poly_clear(reduced_elem);
        fmpz_clear(coeff);

        nf_clear(nf);
    }

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong j;
        nf_t nf;
        nf_elem_t a;
        nmod_poly_t reduced_elem;
        fmpz_t coeff, den;
        ulong mod, d_mod, d_modinv;

        do {
            mod = n_randtest_not_zero(state);
        } while (mod == 1);

        fmpz_init(coeff);
        fmpz_init(den);

        nmod_poly_init(reduced_elem, mod);

        nf_init_randtest(nf, state, 40, 200);

        nf_elem_init(a, nf);

        do {
            nf_elem_randtest(a, state, 200, nf);
            nf_elem_get_den(den, a, nf);
            d_mod = fmpz_fdiv_ui(den, mod);
        } while (n_gcd(d_mod, mod) != 1);

        nf_elem_get_nmod_poly(reduced_elem, a, nf);

        for (j = 0; j < fmpq_poly_degree(nf->pol); j++)
        {
            nf_elem_get_coeff_fmpz(coeff, a, j, nf);
            d_modinv = n_invmod(d_mod, mod);
            result = (nmod_poly_get_coeff_ui(reduced_elem, j) == nmod_mul(fmpz_fdiv_ui(coeff, mod), d_modinv, reduced_elem->mod));
            if (!result)
            {
                printf("FAIL: Reducing element with denominator\n");
                printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
                printf("n = "); flint_printf("%u\n", mod);
                printf("a mod n = "); nmod_poly_print_pretty(reduced_elem, "x"); printf("\n");
                flint_abort();
            }
        }

        fmpz_clear(den);
        fmpz_clear(coeff);
        nf_elem_clear(a, nf);
        nmod_poly_clear(reduced_elem);
        nf_clear(nf);
    }

    TEST_FUNCTION_END(state);
}
