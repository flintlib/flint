/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_elliptic.h"
#include "ca.h"

TEST_FUNCTION_START(ca_elliptic, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_t x, y, ekx, eex, epxy;
        acb_t ay, ekax, eeax, epaxy, aekx, aeex, aepxy;

        ca_ctx_init(ctx);
        ca_init(x, ctx);
        ca_init(y, ctx);
        ca_init(ekx, ctx);
        ca_init(eex, ctx);
        ca_init(epxy, ctx);
        acb_init(ekax);
        acb_init(eeax);
        acb_init(epaxy);
        acb_init(aekx);
        acb_init(aeex);
        acb_init(aepxy);
        acb_init(ay);

        ca_randtest(x, state, 3, 5, ctx);
        if (n_randint(state, 2))
            ca_randtest(y, state, 3, 5, ctx);
        else
            ca_randtest_rational(y, state, 5, ctx);
        if (n_randint(state, 2))
            ca_add(y, y, x, ctx);

        ca_elliptic_k(ekx, x, ctx);
        ca_elliptic_e(eex, x, ctx);
        ca_elliptic_pi(epxy, x, y, ctx);

        ca_get_acb(aekx, ekx, 53, ctx);
        ca_get_acb(aeex, eex, 53, ctx);
        ca_get_acb(aepxy, epxy, 53, ctx);

        ca_get_acb(ekax, x, 53, ctx);
        acb_elliptic_k(ekax, ekax, 53);
        ca_get_acb(eeax, x, 53, ctx);
        acb_elliptic_e(eeax, eeax, 53);
        ca_get_acb(epaxy, x, 53, ctx);
        ca_get_acb(ay, y, 53, ctx);
        acb_elliptic_pi(epaxy, epaxy, ay, 53);

        if (!acb_overlaps(aekx, ekax) || !acb_overlaps(aeex, eeax) ||
            !acb_overlaps(aepxy, epaxy))
        {
            flint_printf("FAIL (Test 1)\n");
            flint_printf("x = "); ca_print(x, ctx); printf("\n\n");
            flint_printf("y = "); ca_print(y, ctx); printf("\n\n");
            flint_printf("ekx = "); ca_print(ekx, ctx); printf("\n\n");
            flint_printf("eex = "); ca_print(eex, ctx); printf("\n\n");
            flint_printf("epxy = "); ca_print(epxy, ctx); printf("\n\n");
            flint_abort();
        }
        
        /* Test the Legendre relation. */
        ca_si_sub(y, 1, x, ctx);
        ca_t eky, eey, p;
        acb_t ekay, eeay, ap, az;
        ca_init(eky, ctx);
        ca_init(eey, ctx);
        ca_init(p, ctx);
        acb_init(ekay);
        acb_init(eeay);
        acb_init(ap);
        acb_init(az);
        acb_zero(az);
        
        ca_pi(p, ctx);
        ca_div_si(p, p, 2, ctx);
        ca_get_acb(ap, p, 53, ctx);
        ca_elliptic_k(eky, y, ctx);
        ca_elliptic_e(eey, y, ctx);
        
        ca_mul(eey, eey, ekx, ctx);
        ca_mul(eex, eex, eky, ctx);
        ca_mul(ekx, ekx, eky, ctx);
        ca_sub(p, p, eey, ctx);
        ca_sub(p, p, eex, ctx);
        ca_add(p, p, ekx, ctx);
        ca_get_acb(ay, p, 53, ctx);
        
        ca_get_acb(ekay, y, 53, ctx);
        acb_elliptic_k(ekay, ekay, 53);
        ca_get_acb(eeay, y, 53, ctx);
        acb_elliptic_e(eeay, eeay, 53);
        
        acb_mul(eeay, eeay, ekax, 53);
        acb_mul(eeax, eeax, ekay, 53);
        acb_mul(ekax, ekax, ekay, 53);
        acb_sub(ap, ap, eeay, 53);
        acb_sub(ap, ap, eeax, 53);
        acb_add(ap, ap, ekax, 53);
        
        if (!acb_overlaps(ap, ay) || !acb_contains(ay, az))
        {
            flint_printf("FAIL (Test 2)\n");
            flint_printf("x = "); ca_print(x, ctx); printf("\n\n");
            flint_printf("y = "); ca_print(y, ctx); printf("\n\n");
            flint_printf("ekx = "); ca_print(ekx, ctx); printf("\n\n");
            flint_printf("eex = "); ca_print(eex, ctx); printf("\n\n");
            flint_printf("eky = "); ca_print(eky, ctx); printf("\n\n");
            flint_printf("eey = "); ca_print(eey, ctx); printf("\n\n");
            flint_printf("ap = %{acb} \n", ap);
            flint_printf("ay = %{acb} \n", ay);
            
            flint_abort();
        }
        
        ca_clear(x, ctx);
        ca_clear(y, ctx);
        ca_clear(ekx, ctx);
        ca_clear(eex, ctx);
        ca_clear(eky, ctx);
        ca_clear(eey, ctx);
        ca_clear(p, ctx);
        ca_clear(epxy, ctx);
        acb_clear(ay);
        acb_clear(ap);
        acb_clear(az);
        acb_clear(ekax);
        acb_clear(eeax);
        acb_clear(ekay);
        acb_clear(eeay);
        acb_clear(epaxy);
        acb_clear(aekx);
        acb_clear(aeex);
        acb_clear(aepxy);
        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
