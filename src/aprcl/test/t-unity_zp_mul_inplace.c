/*
    Copyright (C) 2015 Vladimir Glazachev
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "aprcl.h"

/* unity_zp_mul_inplace and unity_zp_sqr_inplace against the modular
   polynomial arithmetic, for all prime powers p^k used by APRCL and
   larger ones, with reduced and unreduced operands */
TEST_FUNCTION_START(aprcl_unity_zp_mul_inplace, state)
{
    slong i;
    fmpz_t * t;

    t = (fmpz_t *) flint_malloc(sizeof(fmpz_t) * SQUARING_SPACE);
    for (i = 0; i < SQUARING_SPACE; i++)
        fmpz_init(t[i]);

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        ulong p, k, ppow, j;
        fmpz_t n;
        unity_zp f1, f2, g, h;
        fmpz_mod_poly_t cyclo, r;
        int squaring = n_randint(state, 2);

        p = n_randprime(state, 2 + n_randint(state, 5), 0);
        k = (p <= 7) ? 1 + n_randint(state, 3) : 1;
        ppow = n_pow(p, k);

        fmpz_init(n);
        do fmpz_randtest_unsigned(n, state, 1 + n_randint(state, 4000));
        while (fmpz_cmp_ui(n, 1) <= 0);

        unity_zp_init(f1, p, k, n);
        unity_zp_init(f2, p, k, n);
        unity_zp_init(g, p, k, n);
        unity_zp_init(h, p, k, n);

        /* random polynomials, reduced (length <= phi(p^k)) or not; make
           sure to exercise both the full length deg = phi(p^k) (required
           by the straight-line special cases, which read the coefficients
           directly) and, specifically, polynomials truncated by leading
           zero coefficients, which must take the generic path */
        {
            slong deg = (p - 1) * (ppow / p);
            slong glen = 1 + n_randint(state, n_randint(state, 2) ? deg : ppow);
            slong hlen = 1 + n_randint(state, n_randint(state, 2) ? deg : ppow);
            int shape = n_randint(state, 3);

            fmpz_mod_poly_randtest(g->poly, state, glen, g->ctx);
            fmpz_mod_poly_randtest(h->poly, state, hlen, h->ctx);

            if (shape == 1)
            {
                /* full length: leading coefficients forced nonzero */
                fmpz_mod_poly_set_coeff_ui(g->poly, deg - 1, 1 + n_randint(state, 100), g->ctx);
                fmpz_mod_poly_set_coeff_ui(h->poly, deg - 1, 1 + n_randint(state, 100), h->ctx);
            }
            else if (shape == 2 && deg >= 2)
            {
                /* length exactly deg - 1: a leading zero got truncated */
                fmpz_mod_poly_truncate(g->poly, deg - 1, g->ctx);
                fmpz_mod_poly_truncate(h->poly, deg - 1, h->ctx);
                fmpz_mod_poly_set_coeff_ui(g->poly, deg - 2, 1 + n_randint(state, 100), g->ctx);
                fmpz_mod_poly_set_coeff_ui(h->poly, deg - 2, 1 + n_randint(state, 100), h->ctx);
            }
        }

        /* reference: product modulo the cyclotomic polynomial */
        fmpz_mod_poly_init(cyclo, g->ctx);
        fmpz_mod_poly_init(r, g->ctx);
        fmpz_mod_poly_fit_length(cyclo, (p - 1) * (ppow / p) + 1, g->ctx);
        for (j = 0; j < p; j++)
            fmpz_mod_poly_set_coeff_ui(cyclo, j * (ppow / p), 1, g->ctx);
        if (squaring)
            fmpz_mod_poly_mul(r, g->poly, g->poly, g->ctx);
        else
            fmpz_mod_poly_mul(r, g->poly, h->poly, g->ctx);
        fmpz_mod_poly_rem(r, r, cyclo, g->ctx);

        if (squaring)
        {
            unity_zp_sqr_inplace(f1, g, t);
            unity_zp_sqr(f2, g);
        }
        else
        {
            unity_zp_mul_inplace(f1, g, h, t);
            unity_zp_mul(f2, g, h);
        }

        if (!fmpz_mod_poly_equal(f1->poly, r, g->ctx) || !fmpz_mod_poly_equal(f2->poly, r, g->ctx))
        {
            flint_printf("FAIL: p^k = %wu^%wu, squaring = %d\n", p, k, squaring);
            flint_printf("g = "); fmpz_mod_poly_print(g->poly, g->ctx); flint_printf("\n");
            flint_printf("h = "); fmpz_mod_poly_print(h->poly, h->ctx); flint_printf("\n");
            flint_printf("f1 = "); fmpz_mod_poly_print(f1->poly, f1->ctx); flint_printf("\n");
            flint_printf("r = "); fmpz_mod_poly_print(r, g->ctx); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_poly_clear(cyclo, g->ctx);
        fmpz_mod_poly_clear(r, g->ctx);
        unity_zp_clear(f1);
        unity_zp_clear(f2);
        unity_zp_clear(g);
        unity_zp_clear(h);
        fmpz_clear(n);
    }

    for (i = 0; i < SQUARING_SPACE; i++)
        fmpz_clear(t[i]);
    flint_free(t);

    TEST_FUNCTION_END(state);
}
