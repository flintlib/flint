/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(flint_mpn_powmod_preinvn, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        fmpz_t a, d, e, r1, r2;
        nn_ptr dnormed, dinv, as, rs, es;
        mp_size_t n, en;
        flint_bitcnt_t norm, ebits;

        fmpz_init(a);
        fmpz_init(d);
        fmpz_init(e);
        fmpz_init(r1);
        fmpz_init(r2);

        n = 1 + n_randint(state, 12);

        /* any modulus above 1 will do; primality is irrelevant here */
        do {
            fmpz_randbits(d, state, n * FLINT_BITS);
            fmpz_abs(d, d);
        } while (fmpz_size(d) != (slong) n || fmpz_cmp_ui(d, 2) < 0);

        fmpz_randm(a, state, d);

        /* a spread of exponent sizes, so that every window width is used */
        ebits = n_randint(state, 2) ? n_randint(state, 40)
                                    : n_randint(state, 4 * n * FLINT_BITS);
        fmpz_randbits(e, state, ebits);
        fmpz_abs(e, e);

        fmpz_powm(r1, a, e, d);

        en = FLINT_MAX(fmpz_size(e), 1);

        dnormed = flint_malloc((4 * n + en) * sizeof(ulong));
        dinv = dnormed + n;
        as = dinv + n;
        rs = as + n;
        es = rs + n;

        fmpz_get_ui_array(dnormed, n, d);
        norm = flint_clz(dnormed[n - 1]);

        fmpz_get_ui_array(as, n, a);
        fmpz_get_ui_array(es, en, e);

        if (norm)
        {
            mpn_lshift(dnormed, dnormed, n, norm);
            mpn_lshift(as, as, n, norm);
        }

        flint_mpn_preinvn(dinv, dnormed, n);

        flint_mpn_powmod_preinvn(rs, as, es, en, n, dnormed, dinv, norm);

        if (norm)
            mpn_rshift(rs, rs, n, norm);

        fmpz_set_ui_array(r2, rs, n);

        if (!fmpz_equal(r1, r2))
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wd, norm = %wu, ebits = %wu\n",
                    (slong) n, (ulong) norm, (ulong) ebits);
            flint_printf("a = "); fmpz_print(a); flint_printf("\n");
            flint_printf("e = "); fmpz_print(e); flint_printf("\n");
            flint_printf("d = "); fmpz_print(d); flint_printf("\n");
            flint_printf("fmpz_powm = "); fmpz_print(r1); flint_printf("\n");
            flint_printf("flint_mpn = "); fmpz_print(r2); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        flint_free(dnormed);
        fmpz_clear(a);
        fmpz_clear(d);
        fmpz_clear(e);
        fmpz_clear(r1);
        fmpz_clear(r2);
    }

    TEST_FUNCTION_END(state);
}
