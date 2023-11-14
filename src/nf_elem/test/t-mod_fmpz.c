/*
  Copyright (C) 2018 Tommy Hofmann
                2020 Julian Rüth

 ******************************************************************************/

#include "test_helpers.h"
#include "nf_elem.h"

TEST_FUNCTION_START(nf_elem_mod_fmpz, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong j;
        nf_t nf;
        nf_elem_t a, b;
        fmpz_t coeff, mod, reduced_coeff;

        fmpz_init(mod);
        fmpz_randtest_unsigned(mod, state, 2 * FLINT_BITS);
        fmpz_add_ui(mod, mod, 2);

        fmpz_init(coeff);
        fmpz_init(reduced_coeff);

        nf_init_randtest(nf, state, 40, 200);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);

        nf_elem_randtest(a, state, 200, nf);

        nf_elem_mod_fmpz_den(b, a, mod, nf, 0);

        for (j = 0; j < fmpq_poly_degree(nf->pol); j++)
        {
            nf_elem_get_coeff_fmpz(coeff, a, j, nf);
            fmpz_mod(coeff, coeff, mod);
            nf_elem_get_coeff_fmpz(reduced_coeff, b, j, nf);
            result = fmpz_equal(reduced_coeff, coeff);
            if (!result)
            {
                printf("FAIL: Reducing without denominator\n");
                printf("f = "); fmpq_poly_print_pretty(nf->pol, "x"); printf("\n");
                printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
                printf("n = "); fmpz_print(mod); printf("\n");
                printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
                flint_abort();
            }
        }

        nf_elem_smod_fmpz_den(b, a, mod, nf, 0);

        for (j = 0; j < fmpq_poly_degree(nf->pol); j++)
        {
            nf_elem_get_coeff_fmpz(coeff, a, j, nf);
            fmpz_smod(coeff, coeff, mod);
            nf_elem_get_coeff_fmpz(reduced_coeff, b, j, nf);
            result = fmpz_equal(reduced_coeff, coeff);
            if (!result)
            {
                printf("FAIL: Reducing without denominator\n");
                printf("f = "); fmpq_poly_print_pretty(nf->pol, "x"); printf("\n");
                printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
                printf("n = "); fmpz_print(mod); printf("\n");
                printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
                flint_abort();
            }
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        fmpz_clear(coeff);
        fmpz_clear(reduced_coeff);
        fmpz_clear(mod);
        nf_clear(nf);
    }

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong j;
        nf_t nf;
        nf_elem_t a, b, c;
        fmpz_t coeff, mod, den;

        fmpz_init(mod);
        fmpz_randtest_unsigned(mod, state, 2 * FLINT_BITS);
        fmpz_add_ui(mod, mod, 2);

        fmpz_init(coeff);
        fmpz_init(den);

        nf_init_randtest(nf, state, 4, 2);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);

        nf_elem_randtest(a, state, 2, nf);

        nf_elem_mod_fmpz(b, a, mod, nf);

        nf_elem_sub(c, b, a, nf);

        for (j = 0; j < fmpq_poly_degree(nf->pol); j++)
        {
            nf_elem_get_coeff_fmpz(coeff, c, j, nf);
            fmpz_mod(coeff, coeff, mod);
            result = fmpz_is_zero(coeff);
            if (!result || !nf_elem_den_is_one(c, nf))
            {
                printf("FAIL: Reducing without denominator\n");
                printf("f = "); fmpq_poly_print_pretty(nf->pol, "x"); printf("\n");
                printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
                printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
                printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
                printf("n = "); fmpz_print(mod); printf("\n");
                flint_abort();
            }
        }

        nf_elem_smod_fmpz(b, a, mod, nf);

        nf_elem_sub(c, b, a, nf);

        for (j = 0; j < fmpq_poly_degree(nf->pol); j++)
        {
            nf_elem_get_coeff_fmpz(coeff, c, j, nf);
            fmpz_smod(coeff, coeff, mod);
            result = fmpz_is_zero(coeff);
            if (!result || !nf_elem_den_is_one(c, nf))
            {
                printf("FAIL: Reducing without denominator\n");
                printf("f = "); fmpq_poly_print_pretty(nf->pol, "x"); printf("\n");
                printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
                printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
                printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
                printf("n = "); fmpz_print(mod); printf("\n");
                flint_abort();
            }
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);
        fmpz_clear(coeff);
        fmpz_clear(mod);
        fmpz_clear(den);
        nf_clear(nf);
    }

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a, b, c;
        fmpz_t mod, den, gcd;

        fmpz_init(mod);
        fmpz_randtest_unsigned(mod, state, 2 * FLINT_BITS);
        fmpz_add_ui(mod, mod, 2);

        fmpz_init(den);
        fmpz_init(gcd);

        nf_init_randtest(nf, state, 4, 2);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);

        nf_elem_randtest(a, state, 2, nf);

        nf_elem_coprime_den(b, a, mod, nf);

        nf_elem_sub(c, b, a, nf);

        nf_elem_get_den(den, c, nf);
        fmpz_gcd(gcd, den, mod);

        if (!fmpz_is_one(gcd))
        {
            printf("FAIL: Coprime denominators\n");
            printf("f = "); fmpq_poly_print_pretty(nf->pol, "x"); printf("\n");
            printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
            printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
            printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
            printf("n = "); fmpz_print(mod); printf("\n");
            flint_abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);
        fmpz_clear(mod);
        fmpz_clear(den);
        fmpz_clear(gcd);
        nf_clear(nf);
    }

    TEST_FUNCTION_END(state);
}
