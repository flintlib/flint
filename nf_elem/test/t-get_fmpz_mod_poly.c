/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

  Copyright (C) 2018 Tommy Hofmann
                2020 Julian Rüth

 ******************************************************************************/

#include <stdio.h>
#include "nf.h"
#include "nf_elem.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    flint_printf("get_fmpz_mod_poly....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100 * antic_test_multiplier(); i++)
    {
        slong j;
        nf_t nf;
        nf_elem_t a;
        fmpz_mod_poly_t reduced_elem;
        fmpz_t coeff, mod, reduced_coeff;
#if __FLINT_RELEASE >= 20700
        fmpz_mod_ctx_t ctx;
#endif

        fmpz_init(mod);
        fmpz_randtest_unsigned(mod, state, 2 * FLINT_BITS);
        fmpz_add_ui(mod, mod, 2);

        fmpz_init(coeff);
        fmpz_init(reduced_coeff);
#if __FLINT_RELEASE >= 20700
        fmpz_mod_ctx_init(ctx, mod);
        fmpz_mod_poly_init(reduced_elem, ctx);
#else
        fmpz_mod_poly_init(reduced_elem, &mod);
#endif

        nf_init_randtest(nf, state, 40, 200);

        nf_elem_init(a, nf);

        nf_elem_randtest(a, state, 200, nf);

#if __FLINT_RELEASE >= 20700
        nf_elem_get_fmpz_mod_poly_den(reduced_elem, a, nf, 0, ctx);
#else
        nf_elem_get_fmpz_mod_poly_den(reduced_elem, a, nf, 0);
#endif

        for (j = 0; j < fmpq_poly_degree(nf->pol); j++)
        {
            nf_elem_get_coeff_fmpz(coeff, a, j, nf);
            fmpz_mod(coeff, coeff, mod);
#if __FLINT_RELEASE >= 20700
            fmpz_mod_poly_get_coeff_fmpz(reduced_coeff, reduced_elem, j, ctx);
#else
            fmpz_mod_poly_get_coeff_fmpz(reduced_coeff, reduced_elem, j);
#endif
            result = fmpz_equal(reduced_coeff, coeff);
            if (!result)
            {
                printf("FAIL: Reducing without denominator\n");
                printf("f = "); fmpq_poly_print_pretty(nf->pol, "x"); printf("\n");
                printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
                printf("n = "); fmpz_print(mod); printf("\n");
                printf("a mod n = ");
#if __FLINT_RELEASE >= 20700
                fmpz_mod_poly_print_pretty(reduced_elem, "x", ctx);
#else
                fmpz_mod_poly_print_pretty(reduced_elem, "x");
#endif
                printf("\n");
                abort();
            }
        }

        nf_elem_clear(a, nf);
#if __FLINT_RELEASE >= 20700
        fmpz_mod_poly_clear(reduced_elem, ctx);
        fmpz_mod_ctx_clear(ctx);
#else
        fmpz_mod_poly_clear(reduced_elem);
#endif
        fmpz_clear(coeff);
        fmpz_clear(mod);
        nf_clear(nf);
    }

    for (i = 0; i < 100 * antic_test_multiplier(); i++)
    {
        slong j;
        nf_t nf;
        nf_elem_t a;
        fmpz_mod_poly_t reduced_elem;
        fmpz_t coeff, reduced_coeff, den, mod, d_mod, d_modinv;
#if __FLINT_RELEASE >= 20700
        fmpz_mod_ctx_t ctx;
#endif

        fmpz_init(coeff);
        fmpz_init(den);
        fmpz_init(mod);
        fmpz_init(d_mod);
        fmpz_init(d_modinv);
        fmpz_init(reduced_coeff);

        fmpz_init(mod);
        fmpz_randtest_unsigned(mod, state, 2 * FLINT_BITS);
        fmpz_add_ui(mod, mod, 2);

#if __FLINT_RELEASE >= 20700
        fmpz_mod_ctx_init(ctx, mod);
        fmpz_mod_poly_init(reduced_elem, ctx);
#else
        fmpz_mod_poly_init(reduced_elem, &mod);
#endif

        nf_init_randtest(nf, state, 40, 200);

        nf_elem_init(a, nf);

        do {
            nf_elem_randtest(a, state, 200, nf);
            nf_elem_get_den(den, a, nf);
            fmpz_mod(d_mod, den, mod);
            fmpz_gcd(d_mod, d_mod, mod);
        } while (!fmpz_is_one(d_mod));

#if __FLINT_RELEASE >= 20700
        nf_elem_get_fmpz_mod_poly(reduced_elem, a, nf, ctx);
#else
        nf_elem_get_fmpz_mod_poly(reduced_elem, a, nf);
#endif

        for (j = 0; j < fmpq_poly_degree(nf->pol); j++)
        {
            nf_elem_get_coeff_fmpz(coeff, a, j, nf);
            nf_elem_get_den(den, a, nf);
            fmpz_invmod(d_modinv, den, mod);
            fmpz_mul(coeff, coeff, d_modinv);
            fmpz_mod(coeff, coeff, mod);
#if __FLINT_RELEASE >= 20700
            fmpz_mod_poly_get_coeff_fmpz(reduced_coeff, reduced_elem, j, ctx);
#else
            fmpz_mod_poly_get_coeff_fmpz(reduced_coeff, reduced_elem, j);
#endif
            result = (fmpz_equal(coeff, reduced_coeff));
            if (!result)
            {
                printf("FAIL: Reducing element with denominator\n");
                printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
                printf("n = "); fmpz_print(mod); flint_printf("\n");
                printf("a mod n = ");
#if __FLINT_RELEASE >= 20700
                fmpz_mod_poly_print_pretty(reduced_elem, "x", ctx);
#else
                fmpz_mod_poly_print_pretty(reduced_elem, "x");
#endif
                printf("\n");
                abort();
            }
        }

        fmpz_clear(den);
        fmpz_clear(mod);
        fmpz_clear(coeff);
        fmpz_clear(reduced_coeff);
        fmpz_clear(d_mod);
        fmpz_clear(d_modinv);
        nf_elem_clear(a, nf);
#if __FLINT_RELEASE >= 20700
        fmpz_mod_poly_clear(reduced_elem, ctx);
        fmpz_mod_ctx_clear(ctx);
#else
        fmpz_mod_poly_clear(reduced_elem);
#endif
        nf_clear(nf);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
